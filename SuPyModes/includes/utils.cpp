#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include "FinitCoefficient.hpp"

using namespace std;


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



class pBar {
public:
    void update(ScalarType newProgress) {
        currentProgress += newProgress;
        amountOfFiller = (int)((currentProgress / neededProgress)*(ScalarType)pBarLength);
    }
    void print() {
        currUpdateVal %= pBarUpdater.length();
        cout << "\r" //Bring cursor to start of line
            << firstPartOfpBar; //Print out first part of pBar
        for (int a = 0; a < amountOfFiller; a++) { //Print out current progress
            cout << pBarFiller;
        }
        cout << pBarUpdater[currUpdateVal];
        for (int b = 0; b < pBarLength - amountOfFiller; b++) { //Print out spaces
            cout << " ";
        }
        cout << lastPartOfpBar //Print out last part of progress bar
            << " (" << (int)(100*(currentProgress/neededProgress)) << "%)" //This just prints out the percent
            << flush;
        currUpdateVal += 1;
    }
    std::string firstPartOfpBar = "[", //Change these at will (that is why I made them public)
        lastPartOfpBar = "]",
        pBarFiller = "|",
        pBarUpdater = "/-\\|";
private:
    int amountOfFiller,
        pBarLength = 50, //I would recommend NOT changing this
        currUpdateVal = 0; //Do not change
    ScalarType currentProgress = 0, //Do not change
        neededProgress = 100; //I would recommend NOT changing this
};
