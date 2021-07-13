

class EigenSolving{
  public:
    size_t           nMode, MaxIter, Nx, Ny, size;
    uint             DegenerateFactor;
    float            Tolerance, dx, dy, k, lambda, lambdaInit;
    float*           MeshPtr;
    ndarray          Mesh;
    MSparse          Laplacian, EigenMatrix;
    vector<MatrixXf> FullEigenVectors;
    vector<VectorXf> FullEigenValues;

  EigenSolving(ndarray& Mesh,
               size_t   nMode,
               size_t   MaxIter,
               float    Tolerance){
                 this->nMode             = nMode;
                 this->MaxIter           = MaxIter;
                 this->Tolerance         = Tolerance;
                 this->Mesh              = Mesh;
                 this->Nx                = Mesh.request().shape[0];
                 this->Ny                = Mesh.request().shape[1];
                 this->size              = Mesh.request().size;
                 this->dx                = 100./this->Nx;
                 this->dy                = 100./this->Ny;
                 this->MeshPtr           = (float*) Mesh.request().ptr;
               }

   void Setdx(float value){ dx = value;}

   void Setdy(float value){ dy = value;}

   void Setlambda(float value){ lambda = value; k = 2.0 * 3.141592 / lambda;}

   float Getdx(){ return dx;}

   float Getdy(){ return dy;}

   float Getlambda(){ return lambda;}

   ndarray GetMatrix();

   void Points3Laplacian();

   MSparse ComputeMatrix();

   tuple<ndarray, ndarray> PyComputeEigen(float alpha);

   tuple<MatrixXf, VectorXf> ComputeEigen(float alpha);

   void LoopOverITR(ndarray ITRList, float alpha);

   tuple<ndarray, ndarray> GetSlice(uint slice);

   ndarray ComputingOverlap();

   ndarray ComputingCoupling();
};
