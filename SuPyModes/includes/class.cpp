#include "class.hpp"
#include "Extrapolate.hpp"
#define PI 3.1415926535897932384626f
std::complex<ScalarType> J(0.0, 1.0);
typedef std::complex<ScalarType>               ComplexScalarType;


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

    //Vectors.rowwise().reverseInPlace();

    //Values.reverseInPlace();

    return std::make_tuple(Vectors, Values);
    }


void EigenSolving::PringLaplacian(){

  ComputeLaplacian();

  cout<<BaseLaplacian::Laplacian<<endl;
}


void
EigenSolving::ComputeLaplacian(){

  Identity = MSparse(size,size); Identity.setIdentity();

  Points3Laplacian();
}


void
EigenSolving::LoopOverITR(ndarray ITRList, size_t order = 1){


  ITRLength = ITRList.request().size;
  ITRPtr    = (ScalarType*) ITRList.request().ptr;

  this->ITRList    = ITRList;
  this->lambdaInit = lambda;
  this->kInit      = 2.0 * PI / lambda;

  MatrixType EigenVectors, EigenValues;

  ScalarType deltaITR = 0.0;

  FullEigenVectors = vector<MatrixType>(ITRLength);

  FullEigenValues = vector<VectorType>(ITRLength);

  pBar bar;

  ScalarType alpha = -pow( k * ComputeMaxIndex(), 2 );

  for (size_t i=0; i<ITRLength; ++i){

    //bar.update(1/length); bar.print();

    ScalarType lambdaDual = lambdaInit / ITRPtr[i];

    kDual = 2.0 * PI / lambdaDual;

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);

    cout<<"Iteration: "<<i<<endl;

    FullEigenVectors[i] = EigenVectors;

    FullEigenValues[i]  = EigenValues;

    alpha = ExtrapolateNext(order, FullEigenValues, ITRList, i+1);

  }
}




ndarray
EigenSolving::ComputingOverlap(){

  size_t length = FullEigenValues.size()-1,
       iter   = 0;

  VectorType * Overlap = new VectorType(length * nMode * nMode),
             vec0,
             vec1;

  (*Overlap).setZero();

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<nMode; ++i)
          for (size_t j=0; j<nMode; ++j){
              if (j > i){(*Overlap)[iter] = 0.0; ++iter; continue;}
              (*Overlap)[iter] = FullEigenVectors[l].col(i).transpose() * FullEigenVectors[l+1].col(j);
              ++iter; }

  ndarray PyOverlap = Eigen2ndarray( Overlap,
                                     { length, nMode, nMode },
                                     { nMode * nMode * sizeof(ScalarType), nMode * sizeof(ScalarType), sizeof(ScalarType) } );

  return PyOverlap;
}


void
EigenSolving::SortModesFields(){

  size_t length = FullEigenValues.size()-1,
       BestFit,
       iter = 0;

  vector<size_t> Overlap;

  SortedEigenVectors = FullEigenVectors;
  SortedEigenValues = FullEigenValues;

  for (size_t l=0; l<length; ++l){

      Overlap = ComputecOverlaps(SortedEigenVectors[l], FullEigenVectors[l+1]);

        for (size_t j=0; j<nMode; ++j){
            SortedEigenVectors[l+1].col(j) = FullEigenVectors[l+1].col(Overlap[j]);
            SortedEigenValues[l+1][j]        = FullEigenValues[l+1][Overlap[j]];
        }
    }
    FullEigenVectors = SortedEigenVectors;
    FullEigenValues  = SortedEigenValues;
}





void
EigenSolving::SortModesIndex(){

  size_t length = FullEigenValues.size();

  vector<MatrixType> TemporaryVec = FullEigenVectors;

  vector<VectorType> TemporaryVal = FullEigenValues;

  for (size_t l=0; l<length; ++l){
      vector<size_t> sorted = sort_indexes( FullEigenValues[l] );
          for (size_t i=0; i<nMode; ++i)
              TemporaryVec[l].col(i) = FullEigenVectors[l].col(sorted[i]);
      }

  FullEigenVectors = TemporaryVec;

  FullEigenValues  = TemporaryVal;
}


void
EigenSolving::ComputeDegenerateFactor(){
  if (abs(TopSymmetry) == 1)
      DegenerateFactor = 2;
  else
      DegenerateFactor = 1;

  if (abs(BottomSymmetry) == 1)
      DegenerateFactor *= 2;

  if (abs(LeftSymmetry) == 1)
      DegenerateFactor *= 2;

  if (abs(RightSymmetry) == 1)
      DegenerateFactor *= 2;
}

Cndarray
EigenSolving::ComputingCoupling(){

  ComputeDegenerateFactor();

  size_t length = FullEigenValues.size()-1, iter = 0;

  ComplexVectorType * Coupling = new ComplexVectorType(length * nMode * nMode);

  (*Coupling).setZero();

  VectorType vec0, vec1, overlap, temp;

  ScalarType delta, beta0, beta1;

  ComplexScalarType C, I=0.0;

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<nMode; ++i)
          for (size_t j=0; j<nMode; ++j){
              if (j >= i){ (*Coupling)[iter] = 0.0; ++iter; continue; }

              vec0    = FullEigenVectors[l].col(i);
              vec1    = FullEigenVectors[l+1].col(j);

              overlap = vec0.cwiseProduct( vec1 );

              beta0   = FullEigenValues[l][i];
              beta1   = FullEigenValues[l][j];
              delta   = beta0 * beta1;

              temp    = overlap.cwiseProduct( MeshGradient );

              C       = - (ScalarType) 0.5 * J * k*k / sqrt(delta) * abs( 1.0f / (beta0 - beta1) );

              I       = temp.sum();//Trapz(temp, 1.0, Nx, Ny);

              C      *= (ScalarType) DegenerateFactor * I;

              (*Coupling)[iter] = C ;

              ++iter;

            }

  PyCoupling = Eigen2Cndarray( Coupling,
                              { length, nMode, nMode },
                              { nMode * nMode * sizeof(ComplexScalarType), nMode * sizeof(ComplexScalarType), sizeof(ComplexScalarType) } ) ;

  return PyCoupling;

}


ScalarType
EigenSolving::ComputeMaxIndex(){
  MaxIndex = 0.0;
  for (size_t i; i<size; ++i)
     if (MeshPtr[i] > MaxIndex)
         MaxIndex = MeshPtr[i];

 return MaxIndex;
}

Cndarray
EigenSolving::ComputingAdiabatic(){

  ComputeDegenerateFactor();

  size_t iter = 0, length = FullEigenValues.size()-1;

  ComplexVectorType temp;

  VectorType * Adiabatic = new VectorType(length * nMode * nMode);

  (*Adiabatic).setZero();

  VectorType vec0, vec1, overlap;

  ScalarType delta, beta0, beta1;

  ComplexScalarType C, I = 0.0;

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<nMode; ++i)
          for (size_t j=0; j<nMode; ++j){
              if (j >= i){(*Adiabatic)[iter] = 0.0; ++iter; continue;}

              vec0    = FullEigenVectors[l].col(i);
              vec1    = FullEigenVectors[l+1].col(j);

              overlap = vec0.cwiseProduct( vec1 );

              beta0   = FullEigenValues[l][i];
              beta1   = FullEigenValues[l][j];
              delta   = beta0 * beta1;

              temp    = overlap.cwiseProduct( MeshGradient );

              C       = - (ScalarType) 0.5 * J * k*k / sqrt(delta) * abs( 1.0f / (beta0 - beta1) );

              I       = temp.sum();//Trapz(temp, 1.0, Nx, Ny);

              C      *= (ScalarType) DegenerateFactor * I;

              (*Adiabatic)[iter] = abs( delta/C ) ;

              ++iter;
            }

  PyAdiabatic = Eigen2ndarray( Adiabatic,
                              { length, nMode, nMode },
                              { nMode * nMode * sizeof(ScalarType), nMode * sizeof(ScalarType), sizeof(ScalarType) } ) ;

  return PyAdiabatic;

}


tuple<ndarray, ndarray>
EigenSolving::GetSlice(size_t slice){

  MatrixType * Vectors = new MatrixType;
  MatrixType * Values  = new MatrixType;

  (*Vectors) = FullEigenVectors[slice];
  (*Values)  = (-1.0 * FullEigenValues[slice] ).cwiseSqrt() / ITRPtr[slice];

  PyFullEigenVectors = Eigen2ndarray( Vectors,
                                      { nMode, Nx, Ny },
                                      { size * sizeof(ScalarType), Ny * sizeof(ScalarType), sizeof(ScalarType) } ) ;

  PyFullEigenValues = Eigen2ndarray( Values,
                                    { nMode },
                                    { sizeof(ScalarType) } ) ;


  return std::make_tuple( PyFullEigenVectors, PyFullEigenValues );
}




ndarray
EigenSolving::GetFields(){
  size_t slice = 0;
  MatrixType * Vectors = new MatrixType;
  MatrixType * Values  = new MatrixType;

  (*Vectors) = FullEigenVectors[slice];
  (*Values)  = FullEigenValues[slice];

  PyFullEigenVectors = Eigen2ndarray( Vectors,
                                      { nMode, Ny, Nx },
                                      { size * sizeof(ScalarType), Nx * sizeof(ScalarType), sizeof(ScalarType) } ) ;


  return PyFullEigenVectors;

}

ndarray
EigenSolving::GetIndices(){

  ndarray PyIndices(ITRLength * nMode);

  ScalarType  *PyIndicesPtr = (ScalarType*) PyIndices.request().ptr,
              *ITRPtr       = (ScalarType*) ITRList.request().ptr;

  size_t iter = 0;
  cout<<"DDDDDDDDDDDDd"<<kInit<<endl;
  for (size_t l=0; l<ITRLength; ++l)
      for (size_t i=0; i<nMode; ++i){
          PyIndicesPtr[iter] = sqrt(-SortedEigenValues[l][i]) / (ITRPtr[l] * kInit);
          ++iter;}

  PyIndices.resize({ITRLength,nMode});
  return PyIndices;
}






void
BaseLaplacian::LaplacianBoundary(){

  for (size_t j=1; j<size-1; ++j){
      if (j%Ny == 0)
          Laplacian.coeffRef(j-1,j) = 0.0;

      if (j%Ny == Ny-1)
          Laplacian.coeffRef(j+1,j) = 0.0;
  }
}

void
BaseLaplacian::LaplacianXBoundary(){

  int temp = Ny;
  ScalarType dtemp = dy;

  if (LeftSymmetry == 1)
      for (size_t j=0; j<size-1; ++j)
          if (j%temp == 0)
              Laplacian.coeffRef(j,j+1) = 2.0 * SD2A2[0]/pow(dtemp,2);

  if (LeftSymmetry == -1)
      for (size_t j=0; j<size-1; ++j)
          if (j%temp == 0)
              Laplacian.coeffRef(j,j+1) = 0.0;

  if (RightSymmetry == 1)
      for (size_t j=1; j<size; ++j)
          if (j%temp == temp-1)
              Laplacian.coeffRef(j,j-1) = 2.0 * SD2A2[0]/pow(dtemp,2);


  if (RightSymmetry == -1)
      for (size_t j=1; j<size; ++j)
          if (j%temp == temp-1)
              Laplacian.coeffRef(j,j-1) = 0.0;

}

void
BaseLaplacian::LaplacianYBoundary(){

    if (BottomSymmetry == 1)
        for(size_t j=0; j<Ny; ++j)
          Laplacian.coeffRef(j,j+Ny) = 2.0 * SD2A2[0]/pow(dx,2);

    if (BottomSymmetry == -1)
        for(size_t j=0; j<Ny-1; ++j)
            Laplacian.coeffRef(j,j+Ny) = 0.0;

    if (TopSymmetry == 1)
        for(size_t j=size-2*Ny; j<size-Ny; ++j)
            Laplacian.coeffRef(j+Ny,j) = 2.0 * SD2A2[0]/pow(dx,2);

    if (TopSymmetry == -1)
        for(size_t j=size-2*Ny; j<size-Ny; ++j)
            Laplacian.coeffRef(j+Ny,j) = 0.0;
  }


void
BaseLaplacian::Points3Laplacian(){

    Laplacian = MSparse(size,size);

    for(size_t j=0; j<size; ++j)
        Laplacian.insert(j,j) = SD2A2[1]*( 1./pow(dx,2) + 1./pow(dy,2) );

    for(size_t j=0; j<size-1; ++j){
        Laplacian.insert(j+1,j) = SD2A2[0]/pow(dy,2);
        Laplacian.insert(j,j+1) = SD2A2[2]/pow(dy,2);
        }

    for(size_t j=0; j<size-Ny; ++j){
          Laplacian.insert(j+Ny,j) = SD2A2[0]/pow(dx,2);
          Laplacian.insert(j,j+Ny) = SD2A2[2]/pow(dx,2);
          }

    LaplacianBoundary();

    LaplacianXBoundary();

    LaplacianYBoundary();

}





void BaseLaplacian::SetTopSymmetry(int value){
  TopSymmetry = value;
  LaplacianYBoundary();
}

void BaseLaplacian::SetBottomSymmetry(int value){
  BottomSymmetry = value;
  LaplacianYBoundary();
}

void BaseLaplacian::SetLeftSymmetry(int value){
  LeftSymmetry = value;
  LaplacianXBoundary();
}

void BaseLaplacian::SetRightSymmetry(int value){
  RightSymmetry = value;
  LaplacianXBoundary();
}




MSparse
BaseLaplacian::Points5Laplacian(){

  size_t size = Nx*Ny;

  MSparse Laplacian = MSparse(size,size);

  for(size_t j=0; j<size; ++j)
      Laplacian.insert(j,j) = SD2A4[2]*( 1./pow(dx,2) + 1./pow(dx,2) );

  for(size_t j=0; j<size-1; ++j){
      Laplacian.insert(j+1,j) = SD2A4[1]/pow(dx,2);
      Laplacian.insert(j,j+1) = SD2A4[3]/pow(dx,2);
      }


  for(size_t j=0; j<size-2; ++j){
      Laplacian.insert(j+2,j) = SD2A4[0]/pow(dx,2);
      Laplacian.insert(j,j+2) = SD2A4[4]/pow(dx,2);
      }


  for(size_t j=0; j<size-Nx; ++j){
      Laplacian.insert(j+Nx,j) = SD2A4[1]/pow(dy,2);
      Laplacian.insert(j,j+Nx) = SD2A4[3]/pow(dy,2);
      }

  for(size_t j=0; j<size-2*Nx; ++j){
      Laplacian.insert(j+2*Nx,j) = SD2A4[0]/pow(dy,2);
      Laplacian.insert(j,j+2*Nx) = SD2A4[4]/pow(dy,2);
      }

  for(size_t j=2; j<size-2; ++j){
      if (j%Nx == 0)   {Laplacian.coeffRef(j-1,j) = 0.0; Laplacian.coeffRef(j-2,j) = 0.0;}
      if (j%Nx == 1)   {Laplacian.coeffRef(j-2,j) = 0.0;}
      if (j%Nx == Nx-1){Laplacian.coeffRef(j+1,j) = 0.0; Laplacian.coeffRef(j+2,j) = 0.0;}
      if (j%Nx == Nx-2){Laplacian.coeffRef(j+2,j) = 0.0;}

  }


  return Laplacian;
}










// -
