ndarray Eigen2ndarray(MatrixXf *Eigen3Vector, vector<size_t> dimension, vector<size_t> stride){

  ndarray PyVector;

  py::capsule free_when_done(Eigen3Vector->data(), [](void *f) {
     float *foo = reinterpret_cast<float *>(f);
     delete []foo;
   } );

  PyVector = py::array_t<float>( dimension,
                                 stride,
                                 Eigen3Vector->data(),
                                 free_when_done);

   return PyVector;
}

ndarray Eigen2ndarray(VectorXf *Eigen3Vector, vector<size_t> dimension, vector<size_t> stride){

  ndarray PyVector;

  py::capsule free_when_done(Eigen3Vector->data(), [](void *f) {
     float *foo = reinterpret_cast<float *>(f);
     delete []foo;
   } );

  PyVector = py::array_t<float>( dimension,
                                 stride,
                                 Eigen3Vector->data(),
                                 free_when_done);

   return PyVector;
}
