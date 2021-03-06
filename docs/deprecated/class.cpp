#include "class.hpp"
#include "Extrapolate.hpp"
#include "Laplacian.cpp"


using Eigen::internal::BandMatrix;









MSparse
EigenSolving::ComputeMatrix(){

    EigenMatrix = BaseLaplacian::Laplacian;

    size_t iter = 0;

    for(size_t i=0; i<Nx; ++i)
       for(size_t j=0; j<Ny; ++j){
           Identity.coeffRef(iter,iter) = + pow( MeshPtr[iter] * kDual, 2);
           ++iter;
         }

    EigenMatrix += Identity;

    Identity.setIdentity();

    return -1.0*EigenMatrix;
    }


tuple<MatrixType, VectorType>
EigenSolving::ComputeEigen(ScalarType alpha){

    MSparse EigenMatrix = ComputeMatrix();

    SparseGenRealShiftSolve<ScalarType> op(EigenMatrix);

    GenEigsRealShiftSolver<SparseGenRealShiftSolve<ScalarType>> eigs(op, nMode, 2*nMode, alpha);

    eigs.init();

    int nconv = eigs.compute(SortRule::LargestMagn, MaxIter, Tolerance);

    MatrixType Vectors = eigs.eigenvectors().real();

    VectorType Values = eigs.eigenvalues().real();

    return std::make_tuple(Vectors, Values);
    }


void
EigenSolving::ComputeLaplacian(size_t Order){
  Identity = MSparse(size,size); Identity.setIdentity();

  this->Order = Order;

  switch(Order){

  case 2:
      Points3Laplacian(); break;

  case 4:
      Points5Laplacian(); break;
  }

}
























void
EigenSolving::PrepareSuperModes()
{
  for (SuperMode& mode : SuperModes)
      mode.Init(ITRLength, Nx, Ny);
}


void
EigenSolving::PopulateModes(size_t Slice, MatrixType& EigenVectors, VectorType& EigenValues)
{
  for (SuperMode& mode : SuperModes)
  {
    mode.Fields.col(Slice) << EigenVectors.col(mode.ModeNumber);
    mode.Betas[Slice]  = sqrt( - EigenValues[mode.ModeNumber] ) / ITRPtr[Slice];
  }

}


void
EigenSolving::LoopOverITR(ndarray ITRList, size_t order = 1){


  ITRList          = ITRList;
  ITRLength        = ITRList.request().size;
  ITRPtr           = (ScalarType*) ITRList.request().ptr;

  kInit            = 2.0 * PI / lambda;

  MatrixType EigenVectors;
  VectorType EigenValues;




  PrepareSuperModes();


  FullEigenVectors = vector<MatrixType>(ITRLength);
  FullEigenValues = vector<VectorType>(ITRLength);


  ScalarType alpha = -pow( k * ComputeMaxIndex(), 2 );


  for (size_t i=0; i<ITRLength; ++i)
  {

    std::cout<<"Iteration: "<< i <<"   ITR:  "<<ITRPtr[i]<<"\n";

    kDual = kInit * ITRPtr[i] ;

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);

    PopulateModes(i, EigenVectors, EigenValues);

    FullEigenVectors[i] = EigenVectors;

    FullEigenValues[i]  = EigenValues;

    size_t next = i+1, mode=0;

    alpha = ExtrapolateNext(order, FullEigenValues, ITRList, next, mode);
  }
}




ndarray
EigenSolving::GetMode(size_t Mode){

  return SuperModes[Mode].GetFields();
}


ndarray
EigenSolving::GetFullEigen(){

  MatrixType * Vectors = new MatrixType;

  (*Vectors) = Mode0;

  size_t s = Vectors->size();
  ndarray Output = Eigen2ndarray( Vectors, {s} );

  return Output;
}































tuple<VectorType, ScalarType>
EigenSolving::ComputePreSolution(size_t& slice, size_t& mode){

  VectorType X0(size);

  if (slice == 0){
      X0    = VectorType::Ones(size); X0.normalize();
      alpha = pow( k * ComputeMaxIndex(), 2 );
      }
  else{
      X0    = FullEigenVectors[slice-1].col(mode);
      alpha = ExtrapolateNext(ExtrapolOrder, FullEigenValues, ITRList, slice, mode);
      }

  return make_tuple(X0, alpha);
}

void
EigenSolving::PSM(VectorType& X0, size_t& slice, size_t& mode){

  ScalarType ck=0.0, c0=1.0;

  VectorType Y0(size);

  for (size_t iter=0; iter<MaxIter; ++iter){

    GramSchmidt(FullEigenVectors[slice], mode, X0);

    Y0 = Solver.solve(X0);

    ck = Y0.dot(X0);

    Y0.normalize();

    if ( abs( (c0-ck)/ck) < Tolerance && iter!=0 ) break;

    c0 = ck;

    X0 = Y0;
  }

  FullEigenVectors[slice].col(mode)  = Y0;
  FullEigenValues[slice](mode)       = alpha + 1.0/ck;


}






void
EigenSolving::LoopOverITR_(ndarray ITRList, size_t order = 1){

  this->ITRList    = ITRList;
  ITRLength        = ITRList.request().size;
  ITRPtr           = (ScalarType*) ITRList.request().ptr;
  ExtrapolOrder    = order;

  FullEigenVectors = vector<MatrixType>(ITRLength, MatrixType(size, sMode));
  FullEigenValues  = vector<VectorType>(ITRLength, VectorType(sMode));

  kInit            = 2.0 * PI / lambda;

  VectorType X0(size);

  for (size_t slice=0; slice<ITRLength; ++slice){

      std::cout<<"Iteration: "<<slice<<"   ITR:  "<<ITRPtr[slice]<<std::endl;

      kDual = kInit * ITRPtr[slice];

      EigenMatrix = ComputeMatrix();

      for (size_t mode=0; mode<sMode; ++mode){

          tie(X0, alpha) = ComputePreSolution(slice, mode);

          M = -EigenMatrix - (alpha * Identity);

          Solver.compute(M);

          PSM(X0, slice, mode);
        }
    }
}





vector<VectorType>
EigenSolving::ComputeBetas(){
  vector<VectorType> Betas(ITRLength);

  for (size_t i=0; i<ITRLength; ++i)
      Betas[i] = (-FullEigenValues[i]).cwiseSqrt() / ITRPtr[i];

  return Betas;

}



ndarray
EigenSolving::ComputingOverlap(){

  size_t iter   = 0;

  VectorType * Overlap = new VectorType((ITRLength-1) * nMode * nMode),
             vec0,
             vec1;

  (*Overlap).setZero();

  for (size_t l=0; l<ITRLength-1; ++l)
      for (size_t i=0; i<nMode; ++i)
          for (size_t j=0; j<nMode; ++j){
              if (j > i){(*Overlap)[iter] = 0.0; ++iter; continue;}
              (*Overlap)[iter] = FullEigenVectors[l].col(i).transpose() * FullEigenVectors[l+1].col(j);
              ++iter; }

  ndarray PyOverlap = Eigen2ndarray( Overlap, { ITRLength-1, nMode, nMode } );

  return PyOverlap;
}


vector<size_t>
EigenSolving::ComputecOverlaps(MatrixType Matrix0, MatrixType Matrix1, size_t idx){

  size_t nMode = Matrix0.cols();

  ScalarType BestOverlap, Overlap;
  vector<size_t> Indices(sMode);

  for (size_t i=0; i<sMode; ++i){
      BestOverlap = 0;
      for (size_t j=0; j<nMode; ++j){
          Overlap = abs( Matrix0.col(i).transpose() * Matrix1.col(j) );
          if (Overlap > BestOverlap) {Indices[i] = j; BestOverlap = Overlap;}
          }

        if (BestOverlap<0.8)
            std::cout<<"Bad mode correspondence: "<< BestOverlap
                <<"  At ITR: "<< ITRPtr[idx] <<".   You should consider makes more ITR steps"<<std::endl;
      }

  return Indices;

}


void
EigenSolving::SortModes(std::string Type)
{
  std::cout<<"Sorting SuperModes\n";
  if (Type == "Field") SortModesFields();
  else if (Type == "Index") SortModesIndex();
}


void
EigenSolving::SortModesFields(){
  size_t iter = 0;

  vector<size_t> Overlap;

  vector<VectorType> Vec(nMode);
  vector<ScalarType> Val(nMode);

  for (size_t l=0; l<ITRLength-1; ++l){
      Overlap = ComputecOverlaps(FullEigenVectors[l], FullEigenVectors[l+1], l);

        for (size_t j=0; j<nMode; ++j){
            Vec[j] = FullEigenVectors[l+1].col(j);
            Val[j] = FullEigenValues[l+1][j] ;
            }

        for (size_t j=0; j<sMode; ++j){

            FullEigenVectors[l+1].col(j)   = Vec[Overlap[j]];
            FullEigenValues[l+1][j]        = Val[Overlap[j]];
        }
    }
}


void
EigenSolving::SwapMode(SuperMode &Mode0, SuperMode &Mode1)
{
  MatrixType TempField0 = Mode0.Fields,
             TempField1 = Mode1.Fields;

  VectorType TempBeta0 = Mode0.Betas,
             TempBeta1 = Mode1.Betas;

  Mode0.Fields = TempField1;
  Mode1.Fields = TempField0;

  Mode0.Betas = TempBeta1;
  Mode1.Betas = TempBeta0;
}


void
EigenSolving::SortSliceIndex(size_t Slice)
{
  vector<ScalarType> Betas;
  Betas.reserve(nMode);
  vector<SuperMode> TemporaryModes(nMode);

  size_t iter=0;
  for (SuperMode mode : SuperModes)
  {
      Betas.push_back(mode.Betas[Slice]);
      TemporaryModes[iter] = mode;
      ++iter;
  }

  vector<size_t> sorted = sort_indexes( Betas );

  iter = 0;
  for (auto order : sorted)
  {
      SuperModes[iter].CopyOtherSlice(TemporaryModes[order], Slice);
      ++iter;
  }
}


void
EigenSolving::SortModesIndex()
{
  for (size_t l=0; l<ITRLength; ++l)
      SortSliceIndex(l);

}

ndarray
EigenSolving::ComputingCoupling()
{

  std::cout<<"Computing coupling\n";
  VectorType * Coupling = new VectorType(ITRLength, sMode * sMode);

  (*Coupling).setZero();

    size_t iter = 0;
    for (size_t slice=0; slice<ITRLength; ++slice)
        for (size_t i=0; i<sMode; ++i)
            for (size_t j=0; j<sMode; ++j)
                {
                    if (j == i){ (*Coupling)[iter] = 0.0; ++iter; continue; }
                    SuperMode mode0 = SuperModes[i];
                    SuperMode mode1 = SuperModes[j];
                    (*Coupling)[iter] = mode0.ComputeCoupling(mode1, slice, MeshGradient, kInit);
                    ++iter;
                }

    return Eigen2ndarray( Coupling, { ITRLength, sMode, sMode } ) ;
}


ndarray
EigenSolving::ComputingCoupling_(){

  vector<VectorType> Betas = ComputeBetas();

  size_t iter = 0;

  VectorType * Coupling = new VectorType(ITRLength * sMode * sMode);

  (*Coupling).setZero();

  VectorType vec0, vec1, overlap, temp;

  ScalarType delta, beta0, beta1;

  ComplexScalarType C, I;

  for (size_t l=0; l<ITRLength; ++l)
      for (size_t i=0; i<sMode; ++i)
          for (size_t j=0; j<sMode; ++j){
              if (j == i){ (*Coupling)[iter] = 0.0; ++iter; continue; }

              vec0    = FullEigenVectors[l].col(i);
              vec1    = FullEigenVectors[l].col(j);
              beta0   = Betas[l][i];
              beta1   = Betas[l][j];

              overlap = vec0.cwiseProduct( vec1 );

              delta   = beta0 - beta1;

              C       = - (ScalarType) 0.5 * J * kInit*kInit / sqrt(beta0 * beta1) * abs( 1.0f / delta );

              I       = Trapz(overlap.cwiseProduct( MeshGradient ), 1.0, Nx, Ny);

              C      *=  I;

              (*Coupling)[iter] = abs(C);

             ++iter;
            }

  return Eigen2ndarray( Coupling, { ITRLength, sMode, sMode }) ;

}


ScalarType
EigenSolving::ComputeMaxIndex(){
  MaxIndex = 0.0;
  for (size_t i=0; i<size; ++i)
     if (MeshPtr[i] > MaxIndex)
         MaxIndex = MeshPtr[i];

 return MaxIndex;
}



ndarray
EigenSolving::ComputingAdiabatic(){

  vector<VectorType> Betas = ComputeBetas();

  size_t iter = 0;

  VectorType * Adiabatic = new VectorType(ITRLength * sMode * sMode);

  (*Adiabatic).setZero();

  VectorType vec0, vec1, overlap, temp;

  ScalarType delta, beta0, beta1;

  ComplexScalarType C, I;

  for (size_t l=0; l<ITRLength; ++l)
      for (size_t i=0; i<sMode; ++i)
          for (size_t j=0; j<sMode; ++j){
              if (j == i){(*Adiabatic)[iter] = inf; ++iter; continue;}


              vec0    = FullEigenVectors[l].col(i);
              vec1    = FullEigenVectors[l].col(j);
              beta0   = Betas[l][i];
              beta1   = Betas[l][j];

              overlap = vec0.cwiseProduct( vec1 );

              delta   = beta0 - beta1;

              C       = - (ScalarType) 0.5 * J * (kInit*kInit) / sqrt(beta0 * beta1) * abs( 1.0f / delta );

              I       = Trapz(overlap.cwiseProduct( MeshGradient ), 1.0, Nx, Ny);

              C      *=  I;

              (*Adiabatic)[iter] = abs( delta/C ) ;

              ++iter;
            }

  return Eigen2ndarray( Adiabatic, { ITRLength, sMode, sMode } ) ;

}


tuple<ndarray, ndarray>
EigenSolving::GetSlice(size_t slice){

  MatrixType * Vectors = new MatrixType;
  MatrixType * Values  = new MatrixType;

  (*Vectors) = FullEigenVectors[slice];
  (*Values)  = (-1.0 * FullEigenValues[slice] ).cwiseSqrt() / ITRPtr[slice];

  ndarray PyFullEigenVectors = Eigen2ndarray( Vectors, { sMode, Nx, Ny } ) ;

  ndarray PyFullEigenValues = Eigen2ndarray( Values, { sMode } ) ;


  return std::make_tuple( PyFullEigenVectors, PyFullEigenValues );
}





ndarray
EigenSolving::GetFields(size_t slice){

  MatrixType * Vectors = new MatrixType;
  MatrixType * Values  = new MatrixType;

  (*Vectors) = FullEigenVectors[slice];
  (*Values)  = FullEigenValues[slice];

  ndarray PyFullEigenVectors = Eigen2ndarray( Vectors, { sMode, Ny, Nx } ) ;


  return PyFullEigenVectors;

}



ndarray
EigenSolving::GetIndices(){

  ndarray PyIndices(ITRLength * sMode);

  ScalarType  *PyIndicesPtr = (ScalarType*) PyIndices.request().ptr;

  size_t iter = 0;

  for (size_t l=0; l<ITRLength; ++l)
      for (size_t i=0; i<sMode; ++i){
          PyIndicesPtr[iter] = sqrt(abs(FullEigenValues[l][i])) / (ITRPtr[l] * kInit);
          ++iter;}

  PyIndices.resize({ITRLength, sMode});
  return PyIndices;
}


ndarray
EigenSolving::GetBetas(){

  ndarray PyBetas(ITRLength * sMode);

  ScalarType  *PyBetasPtr = (ScalarType*) PyBetas.request().ptr;

  size_t iter = 0;

  for (size_t l=0; l<ITRLength; ++l)
      for (size_t i=0; i<sMode; ++i){
          PyBetasPtr[iter] = sqrt(abs(FullEigenValues[l][i])) / ITRPtr[l];
          ++iter;}

  PyBetas.resize({ITRLength, sMode});

  return PyBetas;
}












// -
