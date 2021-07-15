

class EigenSolving{
  public:
    size_t           nMode, MaxIter, Nx, Ny, size, sMode;
    uint             DegenerateFactor;
    float            Tolerance, dx, dy, k, kInit, lambda, lambdaInit;
    float*           MeshPtr;
    ndarray          Mesh, MeshGradient, ITRList, PyOverlap, PyIndices, PyFullEigenVectors, PyFullEigenValues, PyCoupling;
    MSparse          Laplacian, EigenMatrix;
    vector<MatrixXf> FullEigenVectors, SortedEigenVectors;
    vector<VectorXf> FullEigenValues, SortedEigenValues;

  EigenSolving(ndarray& Mesh,
               //ndarray& MeshGradient,
               size_t   nMode,
               size_t   MaxIter,
               float    Tolerance){
                 this->nMode             = nMode;
                 this->sMode             = nMode;
                 this->MaxIter           = MaxIter;
                 this->Tolerance         = Tolerance;
                 this->Mesh              = Mesh;
                 this->Nx                = Mesh.request().shape[0];
                 this->Ny                = Mesh.request().shape[1];
                 this->size              = Mesh.request().size;
                 this->dx                = 100./this->Nx;
                 this->dy                = 100./this->Ny;
                 this->MeshPtr           = (float*) Mesh.request().ptr;
                 this->DegenerateFactor  = 1.0;
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

   tuple<MatrixXf, VectorXf> ComputeEigen(float alpha);

   void LoopOverITR(ndarray ITRList, float alpha, size_t order);

   tuple<ndarray, ndarray> GetSlice(uint slice);

   ndarray ComputingOverlap();

   Cndarray ComputingCoupling();

   ndarray ComputingIndices();

   void SortModesFields(size_t Mode);

   void SortModesIndex();
};
