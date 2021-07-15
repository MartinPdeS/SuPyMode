#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace std;

vector<size_t> sort_indexes(const VectorXf &v) {

  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}



ndarray
Eigen2ndarray(MatrixXf *Eigen3Vector, vector<size_t> dimension, vector<size_t> stride){

  ndarray PyVector;

  py::capsule free_when_done(Eigen3Vector->data(), [](void *f) {
     float *foo = reinterpret_cast<float *>(f);
     delete []foo; } );

  PyVector = ndarray( dimension, stride, Eigen3Vector->data(), free_when_done );

   return PyVector;
}

ndarray
Eigen2ndarray(VectorXf *Eigen3Vector, vector<size_t> dimension, vector<size_t> stride){

  ndarray PyVector;

  py::capsule free_when_done(Eigen3Vector->data(), [](void *f) {
     float *foo = reinterpret_cast<float *>(f);
     delete []foo; } );

  PyVector = ndarray( dimension, stride, Eigen3Vector->data(), free_when_done );

   return PyVector;
}



Cndarray
Eigen2ndarray(VectorXcf *Eigen3Vector, vector<size_t> dimension, vector<size_t> stride){

  Cndarray PyVector;

  py::capsule free_when_done(Eigen3Vector->data(), [](void *f) {
     fComplex *foo = reinterpret_cast<fComplex *>(f);
     delete []foo;
   } );

  PyVector = Cndarray( dimension, stride, Eigen3Vector->data(), free_when_done);

   return PyVector;
}


float
Trapz(VectorXf& Vector, float dx){
  float sum  = Vector.sum();
  size_t end = Vector.size()-1;
  sum       -= (Vector[0] + Vector[end])/2.0;
  return sum * dx;
}



class pBar {
public:
    void update(double newProgress) {
        currentProgress += newProgress;
        amountOfFiller = (int)((currentProgress / neededProgress)*(double)pBarLength);
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
    double currentProgress = 0, //Do not change
        neededProgress = 100; //I would recommend NOT changing this
};
