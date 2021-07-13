#include "class.hpp"


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

    k = 2.0 * 3.141592 / lambda;

    EigenMatrix = Laplacian;

    MSparse Identity(size,size);

    Identity.setIdentity();

    for(size_t i=0; i<Nx; ++i)
       for(size_t j=0; j<Ny; ++j){
           Identity.coeffRef(i+j*Nx,j*Nx+i) =  2.*(1./pow(dx,2) + 1./pow(dx,2)) - pow( MeshPtr[i+j*Ny]* k, 2);}

    EigenMatrix += Identity;

    return EigenMatrix;
    }


tuple<ndarray, ndarray>
EigenSolving::PyComputeEigen(float alpha){

    MatrixXf * Vectors = new MatrixXf,
             * Values  = new MatrixXf;

    tie(*Vectors, *Values) = ComputeEigen(alpha);

    ndarray EigenVectors = Eigen2ndarray( Vectors, { nMode, Nx, Ny }, { size * sizeof(float), Ny * sizeof(float), sizeof(float) } ) ;

    ndarray EigenValues = Eigen2ndarray( Values, { nMode }, { sizeof(float) } ) ;

    return std::make_tuple(EigenVectors, EigenValues);
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


void
EigenSolving::LoopOverITR(ndarray ITRList, float alpha){

  MatrixXf EigenVectors, EigenValues;

  uint length = ITRList.request().size;

  float* ITRPtr = (float*) ITRList.request().ptr;

  lambdaInit = lambda;

  FullEigenVectors = vector<MatrixXf>(length);

  FullEigenValues = vector<VectorXf>(length);

  for (uint i=0; i<length; ++i){

    lambda = lambdaInit / ITRPtr[i];

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);
    cout<<"Iteration: "<<i<<"  EigenValue:  "<< EigenValues(0) << "  alpha:  " << alpha << "  ITR:  "<<ITRPtr[i]<<endl;

    FullEigenVectors[i] = EigenVectors;

    FullEigenValues[i] = EigenValues;

    alpha = EigenValues(0);// * ITRPtr[i] + EigenValues(nMode-1) * (1-ITRPtr[i]);

  }
}



ndarray
EigenSolving::ComputingOverlap(){

  uint length = FullEigenValues.size()-1;

  VectorXf * Overlap = new VectorXf(length * nMode * nMode);

  VectorXf vec0, vec1;

  uint k = 0;

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<nMode; ++i)
          for (size_t j=0; j<nMode; ++j){
              vec0 = FullEigenVectors[l].col(i);
              vec1 = FullEigenVectors[l+1].col(j);
              (*Overlap)[k] = vec0.transpose() * vec1;
              ++k; }


  return Eigen2ndarray( Overlap,
                        { length, nMode, nMode },
                        { nMode * nMode * sizeof(float), nMode * sizeof(float), sizeof(float) } ) ;

}



ndarray
EigenSolving::ComputingCoupling(){

  uint length = FullEigenValues.size()-1;

  VectorXf * Overlap = new VectorXf(length * nMode * nMode);

  VectorXf vec0, vec1;

  uint k = 0;

  for (size_t l=0; l<length; ++l)
      for (size_t i=0; i<nMode; ++i)
          for (size_t j=0; j<nMode; ++j){
              vec0 = FullEigenVectors[l].col(i);
              vec1 = FullEigenVectors[l+1].col(j);
              (*Overlap)[k] = vec0.transpose() * vec1;
              ++k; }


  return Eigen2ndarray( Overlap,
                        { length, nMode, nMode },
                        { nMode * nMode * sizeof(float), nMode * sizeof(float), sizeof(float) } ) ;

}




tuple<ndarray, ndarray>
EigenSolving::GetSlice(uint slice){

  auto     Vectors = FullEigenVectors[slice];
  auto     Values  = FullEigenValues[slice];

  py::buffer_info InfoVec = py::buffer_info(
                                         Vectors.data(),
                                         Vectors.size() * sizeof(float),
                                         py::format_descriptor<float>::format(),
                                         3,
                                         std::vector<size_t> { nMode, Nx, Ny },
                                         std::vector<size_t> { size * sizeof(float), Ny * sizeof(float), sizeof(float) }
                                         );

   py::buffer_info InfoVal = py::buffer_info(
                                          Values.data(),
                                          Values.size() * sizeof(float),
                                          py::format_descriptor<float>::format(),
                                          1,
                                          std::vector<size_t> { nMode },
                                          std::vector<size_t> { sizeof(float) }
                                          );

    ndarray PyVector = py::array(InfoVec),
            PyValues = py::array(InfoVal);

    return std::make_tuple( PyVector, PyValues );
}
