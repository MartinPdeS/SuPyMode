#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include "FinitCoefficient.hpp"

using namespace std;


VectorType
GetRandomVector(size_t size){
  VectorType Output = VectorType::Random(size);
  Output.normalize();
  return Output;
}

void
GramSchmidt(MatrixType& EigenVectors, size_t& mode)
{

    ScalarType projection;

    for (size_t iter=0; iter<mode; ++iter){

            projection = EigenVectors.col(mode).dot(EigenVectors.col(iter)) / EigenVectors.col(iter).dot(EigenVectors.col(iter));
            EigenVectors.col(mode) = (EigenVectors.col(mode) - projection * EigenVectors.col(iter));
            cout<<"iter: "<<iter<<"   projection:   "<<projection<<endl;
    }
}


void
GramSchmidt(MatrixType& EigenVectors, size_t& mode, VectorType& X0)
{
    if (mode==0)
      return;

    ScalarType projection;

    for (size_t iter=0; iter<mode; ++iter){
            projection = X0.dot(EigenVectors.col(iter)) / EigenVectors.col(iter).dot(EigenVectors.col(iter));
            X0 = X0 - projection * EigenVectors.col(iter);
    }
    X0.normalize();
}

vector<ScalarType>
ComputecOverlaps_(MatrixType EigenVectors0, MatrixType EigenVectors1, size_t iter){

  size_t nMode = EigenVectors0.cols();

  vector<ScalarType> Overlap(nMode);

  for (size_t j=0; j<nMode; ++j)
      Overlap[j] = abs( EigenVectors0.col(iter).transpose() * EigenVectors1.col(j) );

  return Overlap;

}


vector<size_t>
ComputecOverlaps(MatrixType Matrix0, MatrixType Matrix1){

  size_t nMode = Matrix0.cols();

  ScalarType BestOverlap, Overlap;
  vector<size_t> Indices(nMode);

  for (size_t i=0; i<nMode; ++i){
      BestOverlap = 0;
      for (size_t j=0; j<nMode; ++j){
          Overlap = abs( Matrix0.col(i).transpose() * Matrix1.col(j) );
          if (Overlap > BestOverlap) {Indices[i] = j; BestOverlap = Overlap;}
          }
        if (BestOverlap<0.95)
            cout<<"Bad mode correspondence: "<< BestOverlap << "   You should consider makes more ITR steps"<<endl;
      }

  return Indices;

}

vector<size_t>
sort_indexes(const VectorType &v) {

  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


vector<size_t>
sort_indexes(const vector<ScalarType> &v) {

  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}



ndarray
Eigen2ndarray(MatrixType *Eigen3Vector, vector<size_t> dimension, vector<size_t> stride){

  ndarray PyVector;

  py::capsule free_when_done(Eigen3Vector->data(), [](void *f) {
     ScalarType *foo = reinterpret_cast<ScalarType *>(f);
     delete []foo; } );

  PyVector = ndarray( dimension, stride, Eigen3Vector->data(), free_when_done );

   return PyVector;
}

ndarray
Eigen2ndarray(VectorType *Eigen3Vector, vector<size_t> dimension, vector<size_t> stride){

  ndarray PyVector;

  py::capsule free_when_done(Eigen3Vector->data(), [](void *f) {
     ScalarType *foo = reinterpret_cast<ScalarType *>(f);
     delete []foo; } );

  PyVector = ndarray( dimension, stride, Eigen3Vector->data(), free_when_done );

   return PyVector;
}


ndarray
Eigen2ndarray(VectorType &Eigen3Vector, size_t size, vector<size_t> dimension, vector<size_t> stride){


  VectorType * temp = new VectorType(size);

  (*temp) = Eigen3Vector;

  ndarray PyVector;

  py::capsule free_when_done(temp->data(), [](void *f) {
     ScalarType *foo = reinterpret_cast<ScalarType *>(f);
     delete []foo; } );

  PyVector = ndarray( dimension, stride, temp->data(), free_when_done );

   return PyVector;
}



Cndarray
Eigen2Cndarray(ComplexVectorType *Eigen3Vector, vector<size_t> dimension, vector<size_t> stride){

  Cndarray PyVector;

  py::capsule free_when_done(Eigen3Vector->data(), [](void *f) {
     ComplexScalarType *foo = reinterpret_cast<ComplexScalarType *>(f);
     delete []foo;
   } );

  PyVector = Cndarray( dimension, stride, Eigen3Vector->data(), free_when_done);

   return PyVector;
}


ScalarType
Trapz(VectorType& Vector, ScalarType dx, size_t Nx, size_t Ny){

  return Vector.sum();
  ScalarType sum  = 0;
  ScalarType val;

  for (size_t i=0; i<Nx; ++i){
      val = Vector[i];
      if ( i % Nx == 0 || i % Nx == Nx-1 ) { val /= 2.0; };
      if ( i < Ny || i > Ny*(Nx-1) )       { val /= 2.0; };
      sum += val;
      }

  return sum * dx;
}



ScalarType
dydx(ScalarType x, ScalarType y)
{
    return((x - y)/2);
}





//https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
ScalarType
rungeKutta(ScalarType x0, ScalarType y0, ScalarType x, ScalarType h)
{
    // Count number of iterations using step size or
    // step height h
    int n = (int)((x - x0) / h);

    ScalarType k1, k2, k3, k4;

    // Iterate for number of iterations
    ScalarType y = y0;
    for (int i=1; i<=n; i++)
    {
        // Apply Runge Kutta Formulas to find
        // next value of y
        k1 = h*dydx(x0, y);
        k2 = h*dydx(x0 + 0.5*h, y + 0.5*k1);
        k3 = h*dydx(x0 + 0.5*h, y + 0.5*k2);
        k4 = h*dydx(x0 + h, y + k3);

        // Update next value of y
        y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);

        // Update next value of x
        x0 = x0 + h;
    }

    return y;
}
