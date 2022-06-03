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
    mode.Fields.col(Slice)   << EigenVectors.col(mode.ModeNumber);
    mode.Betas[Slice]        = sqrt( - EigenValues[mode.ModeNumber] ) / ITRPtr[Slice];
    mode.EigenValues[Slice]  = EigenValues[mode.ModeNumber];
    mode.Index[Slice]        = sqrt( abs( mode.EigenValues[Slice] ) ) / (ITRPtr[Slice] * kInit);
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


  std::vector<ScalarType> AllFirstEigenValues;
  AllFirstEigenValues.reserve(ITRLength);


  ScalarType alpha = -pow( k * ComputeMaxIndex(), 2 );

  for (size_t i=0; i<ITRLength; ++i)
  {

    std::cout<<"Iteration: "<< i <<"   ITR:  "<<ITRPtr[i]<<"\n";

    kDual = kInit * ITRPtr[i] ;

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);

    PopulateModes(i, EigenVectors, EigenValues);

    AllFirstEigenValues.push_back(EigenValues[0]);

    size_t next = i+1, mode=0;

    alpha = ExtrapolateNext(order, AllFirstEigenValues, ITRList, next);
  }
}




ndarray
EigenSolving::GetMode(size_t Mode){

  return SuperModes[Mode].GetFields();
}
































vector<size_t>
EigenSolving::ComputecOverlaps(MatrixType Matrix0, MatrixType Matrix1, size_t Slice){

  size_t nMode = Matrix0.cols();

  ScalarType BestOverlap, Overlap;
  vector<size_t> Indices(sMode);

  for (size_t i=0; i<sMode; ++i)
    {
      BestOverlap = 0;
      for (size_t j=0; j<nMode; ++j)
          {
              SuperMode Mode0 = SuperModes[i], Mode1 = SuperModes[j];

              Overlap = Mode0.ComputeOverlap(Mode1, Slice);

              if (Overlap > BestOverlap) {Indices[i] = j; BestOverlap = Overlap;}
          }
          if (BestOverlap<0.8)
              std::cout<<"Bad mode correspondence: "<< BestOverlap <<"  At ITR: "<< ITRPtr[Slice] <<". You should consider makes more ITR steps"<<std::endl;
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

  std::cout<<"Computing adiabatic\n";

  VectorType * Adiabatic = new VectorType(ITRLength * sMode * sMode);
  (*Adiabatic).setZero();

  size_t iter = 0;
  for (size_t slice=0; slice<ITRLength; ++slice)
      for (size_t i=0; i<sMode; ++i)
          for (size_t j=0; j<sMode; ++j)
              {
                  if (j == i){ (*Adiabatic)[iter] = 0.0; ++iter; continue; }
                  SuperMode mode0 = SuperModes[i];
                  SuperMode mode1 = SuperModes[j];
                  (*Adiabatic)[iter] = mode0.ComputeAdiabatic(mode1, slice, MeshGradient, kInit);
                  ++iter;
              }

  return Eigen2ndarray( Adiabatic, { ITRLength, sMode, sMode } ) ;

}


tuple<ndarray, ndarray>
EigenSolving::GetSlice(size_t Slice){
  MatrixType OutputFields(size, sMode);
  VectorType OutputBetas(sMode);

  for (size_t mode=0; mode<sMode; ++mode)
  {
    OutputFields.col(mode) = SuperModes[mode].Fields.col(Slice);
    OutputBetas[mode]      = SuperModes[mode].Betas[Slice];
  }


  ndarray FieldsPython = Eigen2ndarray_( OutputFields, { sMode, Nx, Ny } ) ;

  ndarray BetasPython = Eigen2ndarray_( OutputBetas, { sMode } ) ;

  return std::make_tuple( FieldsPython, BetasPython );
}





ndarray
EigenSolving::GetFields(){
  MatrixType Output(ITRLength*size, sMode);

  for (size_t mode=0; mode<sMode; ++mode)
    Output.col(mode)  = SuperModes[mode].Fields;

  return Eigen2ndarray_( Output, { sMode, ITRLength, Nx, Ny } ) ;
}


ndarray
EigenSolving::GetIndices(){
  MatrixType Output(ITRLength, sMode);

  for (size_t mode=0; mode<sMode; ++mode)
    Output.col(mode)  = SuperModes[mode].Index;

  return Eigen2ndarray_( Output, { ITRLength, sMode } ) ;
}


ndarray
EigenSolving::GetBetas(){
  MatrixType Output(ITRLength, sMode);

  for (size_t mode=0; mode<sMode; ++mode)
    Output.col(mode)  = SuperModes[mode].Betas;

  return Eigen2ndarray_( Output, { ITRLength, sMode } ) ;
}












// -
