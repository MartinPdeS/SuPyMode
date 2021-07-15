#include "class.hpp"
#define PI 3.1415926535897932384626f
std::complex<float> J(0.0, 1.0);
typedef std::complex<float>               fComplex;

void
EigenSolving::Points3Laplacian(){

    Laplacian = MSparse(size,size);

    for(uint j=0; j<size; ++j){

        if (j%Nx == 0){
           Laplacian.insert(j+1,j) = -1.0/pow(dx,2);
           continue;
           }

        if (j%Nx == Nx-1){
           Laplacian.insert(j-1,j) = -1.0/pow(dx,2);
           continue;
         }

        Laplacian.insert(j+1,j) = -1.0/pow(dx,2);
        Laplacian.insert(j-1,j) = -1.0/pow(dx,2);
        }

    for(uint j=Nx; j<size; ++j){
        Laplacian.insert(j-Nx,j) = -1.0/pow(dy,2);
        Laplacian.insert(j,j-Nx) = -1.0/pow(dy,2);
        }
    }


ndarray
EigenSolving::GetMatrix(){

    MSparse EigenMatrix = ComputeMatrix();

    MatrixXf temp = MatrixXf(EigenMatrix);

    py::buffer_info Info = py::buffer_info(
                          	               temp.data(),
                          		             temp.size() * sizeof(float),
                          		             py::format_descriptor<float>::format(),
                          		             2,
                          		             std::vector<size_t> {size, size },
                          		             std::vector<size_t> { size * sizeof(float), sizeof(float)}
                                           );

    return py::array(Info);
    }


MSparse
EigenSolving::ComputeMatrix(){

    k = 2.0 * PI / lambda;

    EigenMatrix = Laplacian;

    MSparse Identity(size,size);

    Identity.setIdentity();

    for(size_t i=0; i<Nx; ++i)
       for(size_t j=0; j<Ny; ++j){
           Identity.coeffRef(i+j*Nx,j*Nx+i) =  2.*(1./pow(dx,2) + 1./pow(dx,2)) - pow( MeshPtr[i+j*Ny]* k, 2);}

    EigenMatrix += Identity;

    return EigenMatrix;
    }


tuple<MatrixXf, VectorXf>
EigenSolving::ComputeEigen(float alpha){

    Points3Laplacian();

    MSparse EigenMatrix = ComputeMatrix();

    SparseSymShiftSolve<float> op(EigenMatrix);

    SymEigsShiftSolver<SparseSymShiftSolve<float>> eigs(op, nMode, 2*nMode, alpha);

    eigs.init();

    int nconv = eigs.compute(SortRule::LargestMagn, MaxIter, Tolerance);

    MatrixXf Vectors = eigs.eigenvectors();

    VectorXf Values = eigs.eigenvalues();

    Vectors.rowwise().reverseInPlace();

    Values.reverseInPlace();

    return std::make_tuple(Vectors, Values);
    }



float ExtrapolateNext1(vector<VectorXf>& FullEigenValues, ndarray& ITRList, size_t NextIter){

  return FullEigenValues[NextIter-1][0];
}


float ExtrapolateNext2(vector<VectorXf>& FullEigenValues, ndarray& ITRList, size_t NextIter){


  float * ITRPtr = (float*) ITRList.request().ptr;

  float NextValue;

  float y;
  float y2 = FullEigenValues[NextIter-1][0];
  float y1 = FullEigenValues[NextIter-2][0];

  float x = ITRPtr[NextIter];
  float x2 = ITRPtr[NextIter-1];
  float x1 = ITRPtr[NextIter-2];

  y = (y2 - y1) / (x2 - x1) * (x-x2) + y2;

  return y;
}



float ExtrapolateNext3(vector<VectorXf>& FullEigenValues, ndarray& ITRList, size_t NextIter){


  float * ITRPtr = (float*) ITRList.request().ptr;

  float NextValue;

  float y;
  float y3 = FullEigenValues[NextIter-1][0];
  float y2 = FullEigenValues[NextIter-2][0];
  float y1 = FullEigenValues[NextIter-3][0];

  float x  = ITRPtr[NextIter];
  float x3 = ITRPtr[NextIter-1];
  float x2 = ITRPtr[NextIter-2];
  float x1 = ITRPtr[NextIter-3];

  y = (y2 - y1) / (x2 - x1) * (x-x2) + y2;

  return y;
}



void
EigenSolving::LoopOverITR(ndarray ITRList, float alpha){

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

    lambda = lambdaInit / ITRPtr[i];

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);

    cout<<"Iteration: "<<i<<endl;

    FullEigenVectors[i] = EigenVectors;

    FullEigenValues[i] = EigenValues;

    if (i>1)
        alpha = ExtrapolateNext1(FullEigenValues, ITRList, i+1);

  }
}

ndarray
EigenSolving::ComputingIndices(){

  uint length = SortedEigenValues.size(),
       width  = sMode;

  ndarray PyIndices({length, width});

  float *PyIndicesPtr = (float*) PyIndices.request().ptr,
        *ITRPtr       = (float*) ITRList.request().ptr;

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<width; ++i)
          PyIndicesPtr[l*width + i] = sqrt(-SortedEigenValues[l][i]) / (ITRPtr[l] * kInit);

  return PyIndices;

}


ndarray
EigenSolving::ComputingOverlap(){

  uint length = FullEigenValues.size()-1;

  VectorXf * Overlap = new VectorXf(length * sMode * sMode);

  (*Overlap).setZero();

  VectorXf vec0, vec1;

  uint k = 0;

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<sMode; ++i)
          for (size_t j=0; j<sMode; ++j){
              if (j > i){(*Overlap)[k] = 0.0; ++k; continue;}

              vec0 = FullEigenVectors[l].col(i);
              vec1 = FullEigenVectors[l+1].col(j);
              (*Overlap)[k] = vec0.transpose() * vec1;
              ++k; }


  ndarray PyOverlap = Eigen2ndarray( Overlap,
                                     { length, sMode, sMode },
                                     { sMode * sMode * sizeof(float), sMode * sizeof(float), sizeof(float) } );

  return PyOverlap;
}


void
EigenSolving::SortModesFields(size_t Mode){

  sMode = Mode;

  uint length = FullEigenValues.size()-1,
       BestFit,
       k = 0;

  float  BestOverlap, Overlap;

  VectorXf vec0, vec1;

  SortedEigenVectors = vector<MatrixXf>(length+1);
  SortedEigenValues  = vector<VectorXf>(length+1);


  SortedEigenVectors[0] = MatrixXf(size, sMode);
  SortedEigenValues[0]  = VectorXf(sMode);

  for (size_t i=0; i<sMode; ++i)
      SortedEigenVectors[0].col(i) = FullEigenVectors[0].col(i);

  for (size_t l=0; l<length; ++l){
      SortedEigenVectors[l+1] = MatrixXf(size, sMode);
      SortedEigenValues[l+1] = VectorXf(sMode);

      for (size_t i=0; i<sMode; ++i){
          BestOverlap = 0;
          for (size_t j=0; j<nMode; ++j){

              Overlap = abs( FullEigenVectors[l].col(i).transpose() * FullEigenVectors[l+1].col(j) );

              if ( Overlap > BestOverlap){ BestOverlap = Overlap; BestFit = j; }

              ++k;
            }

            if (BestOverlap < 0.5){cout<<"Bad mode conccordance, you should reconsider the simulations with more steps!   Overlap: "<<BestOverlap;}

            SortedEigenVectors[l+1].col(i) = FullEigenVectors[l+1].col(BestFit);
            SortedEigenValues[l+1][i]      = FullEigenValues[l+1][BestFit];
        }
    }
}


void
EigenSolving::SortModesIndex(){

  uint length = FullEigenValues.size();

  vector<MatrixXf> TemporaryVec = FullEigenVectors;

  vector<VectorXf> TemporaryVal = FullEigenValues;

  for (size_t l=0; l<length; ++l){
      vector<size_t> sorted = sort_indexes( FullEigenValues[l] );
          for (size_t i=0; i<sMode; ++i)
              TemporaryVec[l].col(i) = FullEigenVectors[l].col(sorted[i]);
      }

  FullEigenVectors = TemporaryVec;

  FullEigenValues  = TemporaryVal;
}



Cndarray
EigenSolving::ComputingCoupling(){

  uint length = FullEigenValues.size()-1;

  VectorXcf * Coupling = new VectorXcf(length * sMode * sMode);

  (*Coupling).setZero();

  VectorXf vec0, vec1;

  uint k = 0;

  VectorXf overlap;

  float delta, beta0, beta1, I=0.0;

  fComplex C;

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<sMode; ++i)
          for (size_t j=0; j<sMode; ++j){
              if (j > i){(*Coupling)[k] = 0.0; ++k; continue;}

              vec0    = FullEigenVectors[l].col(i);
              vec1    = FullEigenVectors[l+1].col(j);

              overlap = vec0.cwiseProduct( vec1 );

              beta0   = FullEigenValues[l][i];
              beta1   = FullEigenValues[l][j];
              delta   = beta0 * beta1;

              C       = -0.5 * k*k / sqrt(delta) * abs( 1.0 / (beta0 - beta1) );

              I       = Trapz(overlap, 1.0);

              C      *= DegenerateFactor * I;

              (*Coupling)[k] = C * J;

              ++k;

            }

  PyCoupling = Eigen2ndarray( Coupling,
                              { length, sMode, sMode },
                              { sMode * sMode * sizeof(fComplex), sMode * sizeof(fComplex), sizeof(fComplex) } ) ;

  return PyCoupling;

}



tuple<ndarray, ndarray>
EigenSolving::GetSlice(uint slice){

  MatrixXf * Vectors = new MatrixXf;
  MatrixXf * Values  = new MatrixXf;

  (*Vectors) = SortedEigenVectors[slice];
  (*Values)  = SortedEigenValues[slice];

  PyFullEigenVectors = Eigen2ndarray( Vectors,
                                      { sMode, Nx, Ny },
                                      { size * sizeof(float), Ny * sizeof(float), sizeof(float) } ) ;

  PyFullEigenValues = Eigen2ndarray( Values,
                                    { sMode },
                                    { sizeof(float) } ) ;


  return std::make_tuple( PyFullEigenVectors, PyFullEigenValues );
}


















// -
