#include "class.hpp"
#include "Extrapolate.hpp"
#define PI 3.1415926535897932384626f
std::complex<float> J(0.0, 1.0);
typedef std::complex<float>               fComplex;


MSparse
EigenSolving::ComputeMatrix(){

    EigenMatrix = BaseLaplacian::Laplacian;
    uint iter = 0;

    for(size_t i=0; i<Nx; ++i)
       for(size_t j=0; j<Ny; ++j){
           Identity.coeffRef(iter,iter) = + pow( MeshPtr[iter] * kDual, 2);
           ++iter;
         }

    EigenMatrix += Identity;

    return -1.0*EigenMatrix;
    }


tuple<MatrixXf, VectorXf>
EigenSolving::ComputeEigen(float alpha){

    MSparse EigenMatrix = ComputeMatrix();

    SparseGenRealShiftSolve<float> op(EigenMatrix);

    Spectra::GenEigsRealShiftSolver<SparseGenRealShiftSolve<float>> eigs(op, nMode, 2*nMode, alpha);

    eigs.init();

    int nconv = eigs.compute(SortRule::LargestMagn, MaxIter, Tolerance);

    MatrixXf Vectors = eigs.eigenvectors().real();

    VectorXf Values = eigs.eigenvalues().real();

    Vectors.rowwise().reverseInPlace();

    Values.reverseInPlace();

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
EigenSolving::LoopOverITR(ndarray ITRList, float alpha, size_t order = 1){


  uint length = ITRList.request().size;
  float* ITRPtr = (float*) ITRList.request().ptr;

  this->ITRList    = ITRList;
  this->lambdaInit = lambda;
  this->kInit      = 2.0 * PI / lambda;

  MatrixXf EigenVectors, EigenValues;

  float deltaITR = 0.0;


  FullEigenVectors = vector<MatrixXf>(length);

  FullEigenValues = vector<VectorXf>(length);

  pBar bar;

  for (uint i=0; i<length; ++i){

    //bar.update(1/length); bar.print();

    float lambdaDual = lambdaInit / ITRPtr[i];

    kDual = 2.0 * PI / lambdaDual;

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);

    cout<<"Iteration: "<<i<<endl;

    FullEigenVectors[i] = EigenVectors;

    FullEigenValues[i]  = EigenValues;

    if (i>order)
        alpha = ExtrapolateNext(order, FullEigenValues, ITRList, i+1);

  }
}

ndarray
EigenSolving::ComputingIndices(){

  uint    length       = SortedEigenValues.size(),
          width        = nMode;

  ndarray PyIndices({length, width});

  float  *PyIndicesPtr = (float*) PyIndices.request().ptr,
         *ITRPtr       = (float*) ITRList.request().ptr;

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<width; ++i)
          PyIndicesPtr[l*width + i] = sqrt(-SortedEigenValues[l][i]) / (ITRPtr[l] * kInit);

  return PyIndices;

}


ndarray
EigenSolving::ComputingOverlap(){

  uint length = FullEigenValues.size()-1,
       iter   = 0;

  VectorXf * Overlap = new VectorXf(length * nMode * nMode),
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
                                     { nMode * nMode * sizeof(float), nMode * sizeof(float), sizeof(float) } );

  return PyOverlap;
}


void
EigenSolving::SortModesFields(){

  uint length = FullEigenValues.size()-1,
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

  uint length = FullEigenValues.size();

  vector<MatrixXf> TemporaryVec = FullEigenVectors;

  vector<VectorXf> TemporaryVal = FullEigenValues;

  for (size_t l=0; l<length; ++l){
      vector<size_t> sorted = sort_indexes( FullEigenValues[l] );
          for (size_t i=0; i<nMode; ++i)
              TemporaryVec[l].col(i) = FullEigenVectors[l].col(sorted[i]);
      }

  FullEigenVectors = TemporaryVec;

  FullEigenValues  = TemporaryVal;
}



Cndarray
EigenSolving::ComputingCoupling(){

  uint length = FullEigenValues.size()-1, iter = 0;

  VectorXcf * Coupling = new VectorXcf(length * nMode * nMode);

  (*Coupling).setZero();

  VectorXf vec0, vec1, overlap, temp;

  float delta, beta0, beta1;

  fComplex C, I=0.0;

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

              C       = -0.5f * J * k*k / sqrt(delta) * abs( 1.0f / (beta0 - beta1) );

              I       = temp.sum();//Trapz(temp, 1.0, Nx, Ny);

              C      *= (float) DegenerateFactor * I;

              (*Coupling)[iter] = C ;

              ++iter;

            }

  PyCoupling = Eigen2Cndarray( Coupling,
                              { length, nMode, nMode },
                              { nMode * nMode * sizeof(fComplex), nMode * sizeof(fComplex), sizeof(fComplex) } ) ;

  return PyCoupling;

}


Cndarray
EigenSolving::ComputingAdiabatic(){

  uint iter = 0, length = FullEigenValues.size()-1;

  VectorXcf temp;

  VectorXf * Adiabatic = new VectorXf(length * nMode * nMode);

  (*Adiabatic).setZero();

  VectorXf vec0, vec1, overlap;

  float delta, beta0, beta1;

  fComplex C, I = 0.0;

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

              C       = -0.5f * J * k*k / sqrt(delta) * abs( 1.0f / (beta0 - beta1) );

              I       = temp.sum();//Trapz(temp, 1.0, Nx, Ny);

              C      *= (float) DegenerateFactor * I;

              (*Adiabatic)[iter] = abs( delta/C ) ;

              ++iter;
            }

  PyAdiabatic = Eigen2ndarray( Adiabatic,
                              { length, nMode, nMode },
                              { nMode * nMode * sizeof(float), nMode * sizeof(float), sizeof(float) } ) ;

  return PyAdiabatic;

}


tuple<ndarray, ndarray>
EigenSolving::GetSlice(uint slice){

  MatrixXf * Vectors = new MatrixXf;
  MatrixXf * Values  = new MatrixXf;

  (*Vectors) = FullEigenVectors[slice];
  (*Values)  = FullEigenValues[slice];

  PyFullEigenVectors = Eigen2ndarray( Vectors,
                                      { nMode, Ny, Nx },
                                      { size * sizeof(float), Nx * sizeof(float), sizeof(float) } ) ;

  PyFullEigenValues = Eigen2ndarray( Values,
                                    { nMode },
                                    { sizeof(float) } ) ;


  return std::make_tuple( PyFullEigenVectors, PyFullEigenValues );
}















void
BaseLaplacian::LaplacianBoundary(){

  for (uint j=1; j<size-1; ++j){
      if (j%Ny == 0)
          Laplacian.coeffRef(j-1,j) = 0.0;

      if (j%Ny == Ny-1)
          Laplacian.coeffRef(j+1,j) = 0.0;
  }
}

void
BaseLaplacian::LaplacianXBoundary(){

  int temp = Ny;
  float dtemp = dy;

  if (LeftSymmetry == 1)
      for (uint j=0; j<size-1; ++j)
          if (j%temp == 0)
              Laplacian.coeffRef(j,j+1) = 2.0 * SD2A2[0]/pow(dtemp,2);

  if (LeftSymmetry == -1)
      for (uint j=0; j<size-1; ++j)
          if (j%temp == 0)
              Laplacian.coeffRef(j,j+1) = 0.0;

  if (RightSymmetry == 1)
      for (uint j=1; j<size; ++j)
          if (j%temp == temp-1)
              Laplacian.coeffRef(j,j-1) = 2.0 * SD2A2[0]/pow(dtemp,2);


  if (RightSymmetry == -1)
      for (uint j=1; j<size; ++j)
          if (j%temp == temp-1)
              Laplacian.coeffRef(j,j-1) = 0.0;

}

void
BaseLaplacian::LaplacianYBoundary(){

    if (BottomSymmetry == 1)
        for(uint j=0; j<Ny; ++j)
          Laplacian.coeffRef(j,j+Ny) = 2.0 * SD2A2[0]/pow(dx,2);

    if (BottomSymmetry == -1)
        for(uint j=0; j<Ny-1; ++j)
            Laplacian.coeffRef(j,j+Ny) = 0.0;

    if (TopSymmetry == 1)
        for(uint j=size-2*Ny; j<size-Ny; ++j)
            Laplacian.coeffRef(j+Ny,j) = 2.0 * SD2A2[0]/pow(dx,2);

    if (TopSymmetry == -1)
        for(uint j=size-2*Ny; j<size-Ny; ++j)
            Laplacian.coeffRef(j+Ny,j) = 0.0;
  }


void
BaseLaplacian::Points3Laplacian(){

    Laplacian = MSparse(size,size);

    for(uint j=0; j<size; ++j)
        Laplacian.insert(j,j) = SD2A2[1]*( 1./pow(dx,2) + 1./pow(dy,2) );

    for(uint j=0; j<size-1; ++j){
        Laplacian.insert(j+1,j) = SD2A2[0]/pow(dy,2);
        Laplacian.insert(j,j+1) = SD2A2[2]/pow(dy,2);
        }

    for(uint j=0; j<size-Ny; ++j){
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

  for(uint j=0; j<size; ++j)
      Laplacian.insert(j,j) = SD2A4[2]*( 1./pow(dx,2) + 1./pow(dx,2) );

  for(uint j=0; j<size-1; ++j){
      Laplacian.insert(j+1,j) = SD2A4[1]/pow(dx,2);
      Laplacian.insert(j,j+1) = SD2A4[3]/pow(dx,2);
      }


  for(uint j=0; j<size-2; ++j){
      Laplacian.insert(j+2,j) = SD2A4[0]/pow(dx,2);
      Laplacian.insert(j,j+2) = SD2A4[4]/pow(dx,2);
      }


  for(uint j=0; j<size-Nx; ++j){
      Laplacian.insert(j+Nx,j) = SD2A4[1]/pow(dy,2);
      Laplacian.insert(j,j+Nx) = SD2A4[3]/pow(dy,2);
      }

  for(uint j=0; j<size-2*Nx; ++j){
      Laplacian.insert(j+2*Nx,j) = SD2A4[0]/pow(dy,2);
      Laplacian.insert(j,j+2*Nx) = SD2A4[4]/pow(dy,2);
      }

  for(uint j=2; j<size-2; ++j){
      if (j%Nx == 0)   {Laplacian.coeffRef(j-1,j) = 0.0; Laplacian.coeffRef(j-2,j) = 0.0;}
      if (j%Nx == 1)   {Laplacian.coeffRef(j-2,j) = 0.0;}
      if (j%Nx == Nx-1){Laplacian.coeffRef(j+1,j) = 0.0; Laplacian.coeffRef(j+2,j) = 0.0;}
      if (j%Nx == Nx-2){Laplacian.coeffRef(j+2,j) = 0.0;}

  }


  return Laplacian;
}










// -
