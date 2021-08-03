#include "class.hpp"
#include "Extrapolate.hpp"
#include "Laplacian.cpp"

std::complex<ScalarType> J(0.0, 1.0);
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

    if (Debug)
        cout<<EigenMatrix<<endl;

    return -1.0*EigenMatrix;
    }


tuple<MatrixType, VectorType>
EigenSolving::ComputeEigen(ScalarType alpha){

    MSparse EigenMatrix = ComputeMatrix();

    SparseGenRealShiftSolve<ScalarType> op(EigenMatrix);

    Spectra::GenEigsRealShiftSolver<SparseGenRealShiftSolve<ScalarType>> eigs(op, nMode, 2*nMode, alpha);

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
EigenSolving::LoopOverITR(ndarray ITRList, size_t order = 1){


  ITRLength = ITRList.request().size;
  ITRPtr    = (ScalarType*) ITRList.request().ptr;

  this->ITRList    = ITRList;
  this->lambdaInit = lambda;
  this->kInit      = 2.0 * PI / lambda;

  MatrixType EigenVectors, EigenValues;

  FullEigenVectors = vector<MatrixType>(ITRLength);

  FullEigenValues = vector<VectorType>(ITRLength);

  ScalarType alpha = -pow( k * ComputeMaxIndex(), 2 );

  for (size_t i=0; i<ITRLength; ++i){

    ScalarType lambdaDual = lambdaInit / ITRPtr[i];

    kDual = 2.0 * PI / lambdaDual;

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);

    cout<<"Iteration: "<<i<<"   ITR:  "<<ITRPtr[i]<<endl;

    FullEigenVectors[i] = EigenVectors;

    FullEigenValues[i]  = EigenValues;

    alpha = ExtrapolateNext(order, FullEigenValues, ITRList, i+1);
  }
}




tuple<VectorType, ScalarType>
EigenSolving::PSM(ConjugateGradient<MSparse>& solver, VectorType& X0){

  VectorType Y0(size);

  ScalarType ck, c0;

  for (size_t iter=0; iter<MaxIter; ++iter){

    Y0 = solver.solve(X0);

    ck = Y0.dot(X0)/X0.dot(X0);

    Y0.normalize();

    X0 = Y0;

    if (abs(ck-c0) < Tolerance)
        break;

    c0 = ck;

  }

  return make_tuple(X0,ck) ;

}

ndarray
EigenSolving::LoopOverITR_(ndarray ITRList, size_t order = 1, ScalarType lol=0.0){

  ndarray output;

  VectorType * temp = new VectorType(size);

  kDual = 2.0 * PI / lambda;

  MSparse EigenMatrix = ComputeMatrix();

  VectorType EigenVector;

  ScalarType alpha = pow( k * ComputeMaxIndex(), 2 ), EigenValue;

  Identity.setIdentity();

  MSparse M = -EigenMatrix - (alpha * Identity);

  ConjugateGradient<MSparse> solver(M);

  VectorType X0(size); X0.setOnes(); X0.normalize();

  tie(EigenVector, EigenValue) = PSM(solver, X0);

  EigenValue = 1.0 / EigenValue + alpha;

  (*temp) = EigenVector;

  ndarray PyOutput = Eigen2ndarray( temp,
                                    { Nx, Ny },
                                    { Nx * sizeof(ScalarType), sizeof(ScalarType) } );

  cout<<EigenValue<<endl;
  return PyOutput;
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

  ndarray PyOverlap = Eigen2ndarray( Overlap,
                                     { ITRLength-1, nMode, nMode },
                                     { nMode * nMode * sizeof(ScalarType), nMode * sizeof(ScalarType), sizeof(ScalarType) } );

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
            cout<<"Bad mode correspondence: "<< BestOverlap
                <<"  At ITR: "<< ITRPtr[idx] <<".   You should consider makes more ITR steps"<<endl;
      }

  return Indices;

}

void
EigenSolving::SortModesFields(){
  size_t BestFit, iter = 0;

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
EigenSolving::SortModesIndex(){

  SortedEigenVectors = FullEigenVectors;

  SortedEigenValues = FullEigenValues;

  for (size_t l=0; l<ITRLength; ++l){
      vector<size_t> sorted = sort_indexes( FullEigenValues[l] );
          for (size_t i=0; i<nMode; ++i){
              SortedEigenVectors[l].col(i)   = FullEigenVectors[l].col(sorted[i]);
              SortedEigenValues[l][i]        = FullEigenValues[l][ sorted[i] ];
        }
    }

  FullEigenVectors = SortedEigenVectors;

  FullEigenValues  = SortedEigenValues;
}


void
EigenSolving::ComputeDegenerateFactor(){

  DegenerateFactor = 1;

  if (abs(TopSymmetry) == 1)
      DegenerateFactor = 2;

  if (abs(LeftSymmetry) == 1)
      DegenerateFactor *= 2;

  if (abs(RightSymmetry) == 1)
      DegenerateFactor *= 2;

  if (abs(BottomSymmetry) == 1)
      DegenerateFactor *= 2;

}


ndarray
EigenSolving::ComputingCoupling(){

  vector<VectorType> Betas = ComputeBetas();

  ComputeDegenerateFactor();

  size_t iter = 0;

  VectorType * Coupling = new VectorType((ITRLength-1) * sMode * sMode);

  (*Coupling).setZero();

  VectorType vec0, vec1, overlap, temp;

  ScalarType delta, beta0, beta1;

  ComplexScalarType C, I;

  for (size_t l=1; l<ITRLength; ++l)
      for (size_t i=0; i<sMode; ++i)
          for (size_t j=0; j<sMode; ++j){
              if (j == i){ (*Coupling)[iter] = 0.0; ++iter; continue; }

              vec0    = FullEigenVectors[l].col(i);
              vec1    = FullEigenVectors[l].col(j);
              beta0   = Betas[l][i];
              beta1   = Betas[l][j];

              overlap = vec0.cwiseProduct( vec1 );

              //if (overlap.cwiseAbs().sum() > 0.9){ (*Coupling)[iter] = 0.0; ++iter; continue; }

              delta   = beta0 - beta1;

              temp    = overlap.cwiseProduct( MeshGradient );

              C       = - (ScalarType) 0.5 * J * kInit*kInit / sqrt(beta0 * beta1) * abs( 1.0f / delta );

              I       = Trapz(temp, 1.0, Nx, Ny);

              //C      *= (ComplexScalarType) DegenerateFactor * I;

              C      *=  I;


              (*Coupling)[iter] = abs(C);

             ++iter;
            }

  return Eigen2ndarray( Coupling,
                        { ITRLength-1, sMode, sMode },
                        { sMode * sMode * sizeof(ScalarType), sMode * sizeof(ScalarType), sizeof(ScalarType) } ) ;

}


ScalarType
EigenSolving::ComputeMaxIndex(){
  MaxIndex = 0.0;
  for (size_t i; i<size; ++i)
     if (MeshPtr[i] > MaxIndex)
         MaxIndex = MeshPtr[i];

 return MaxIndex;
}



ndarray
EigenSolving::ComputingAdiabatic(){

  ITRPtr    = (ScalarType*) ITRList.request().ptr;

  vector<VectorType> Betas = ComputeBetas();

  ComputeDegenerateFactor();

  size_t iter = 0;

  VectorType temp;

  VectorType * Adiabatic = new VectorType((ITRLength-1) * sMode * sMode);

  (*Adiabatic).setZero();

  VectorType vec0, vec1, overlap;

  ScalarType delta, beta0, beta1;

  ComplexScalarType C, I = 0.0;

  for (size_t l=1; l<ITRLength; ++l)
      for (size_t i=0; i<sMode; ++i)
          for (size_t j=0; j<sMode; ++j){
              if (j == i){(*Adiabatic)[iter] = inf; ++iter; continue;}


              vec0    = FullEigenVectors[l].col(i);
              vec1    = FullEigenVectors[l].col(j);
              beta0   = Betas[l][i];
              beta1   = Betas[l][j];
              delta   = beta0 - beta1;

              overlap = vec0.cwiseProduct( vec1 );

              //if (overlap.cwiseAbs().sum() > 0.9){ (*Adiabatic)[iter] = inf; ++iter; continue; }

              temp    = overlap.cwiseProduct( MeshGradient );

              C       = - (ScalarType) 0.5 * J * (kInit*kInit) / sqrt(beta0 * beta1) * abs( 1.0f / delta );

              I       = Trapz(temp, 1.0, Nx, Ny);

              //C      *= (ComplexScalarType) DegenerateFactor * I;

              C      *=  I;

              (*Adiabatic)[iter] = abs( delta/C ) ;

              ++iter;
            }

  return Eigen2ndarray( Adiabatic,
                        { ITRLength-1, sMode, sMode },
                        { sMode * sMode * sizeof(ScalarType), sMode * sizeof(ScalarType), sizeof(ScalarType) } ) ;

}


tuple<ndarray, ndarray>
EigenSolving::GetSlice(size_t slice){

  MatrixType * Vectors = new MatrixType;
  MatrixType * Values  = new MatrixType;

  (*Vectors) = this->FullEigenVectors[slice];
  (*Values)  = (-1.0 * this->FullEigenValues[slice] ).cwiseSqrt() / ITRPtr[slice];

  ndarray PyFullEigenVectors = Eigen2ndarray( Vectors,
                                      { sMode, Nx, Ny },
                                      { size * sizeof(ScalarType), Ny * sizeof(ScalarType), sizeof(ScalarType) } ) ;

  ndarray PyFullEigenValues = Eigen2ndarray( Values,
                                             { sMode },
                                             { sizeof(ScalarType) } ) ;


  return std::make_tuple( PyFullEigenVectors, PyFullEigenValues );
}




ndarray
EigenSolving::GetFields(size_t slice){

  MatrixType * Vectors = new MatrixType;
  MatrixType * Values  = new MatrixType;

  (*Vectors) = FullEigenVectors[slice];
  (*Values)  = FullEigenValues[slice];

  ndarray PyFullEigenVectors = Eigen2ndarray( Vectors,
                                      { sMode, Ny, Nx },
                                      { size * sizeof(ScalarType), Nx * sizeof(ScalarType), sizeof(ScalarType) } ) ;


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
