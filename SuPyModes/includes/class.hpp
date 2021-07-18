class BaseLaplacian{
  public:
    int      TopSymmetry, BottomSymmetry, LeftSymmetry, RightSymmetry;
    ndarray  Mesh;
    size_t   Nx, Ny, size;
    float    dx, dy;
    bool     debug;
    MSparse  Laplacian;

    BaseLaplacian(ndarray&  Mesh){
      this->Nx                = Mesh.request().shape[0];
      this->Ny                = Mesh.request().shape[1];
      this->size              = Mesh.request().size;
      this->TopSymmetry      = 0;
      this->BottomSymmetry   = 0;
      this->LeftSymmetry     = 0;
      this->RightSymmetry    = 0;
      this->debug             = false;
    }

    void SetLeftSymmetry(int value);
    int GetLeftSymmetry(){ return LeftSymmetry; }

    void SetRightSymmetry(int value);
    int GetRightSymmetry(){ return RightSymmetry; }

    void SetTopSymmetry(int value);
    int GetTopSymmetry(){ return TopSymmetry; }

    void SetBottomSymmetry(int value);
    int GetBottomSymmetry(){ return BottomSymmetry; }

    void Points3Laplacian();

    void LaplacianBoundary();

    void LaplacianXBoundary();

    void LaplacianYBoundary();

    MSparse Points5Laplacian();

    void Setdx(float value){ dx = value; }

    void Setdy(float value){ dy = value; }
};


class EigenSolving : public BaseLaplacian{
  public:
    bool             debug;
    size_t           nMode, MaxIter, Nx, Ny, size, sMode;
    uint             DegenerateFactor;
    float            Tolerance, dx, dy, k, kInit, kDual, lambda, lambdaInit;
    float*           MeshPtr;
    ndarray          Mesh, ITRList, PyOverlap, PyIndices, PyFullEigenVectors, PyFullEigenValues, PyAdiabatic;
    Cndarray         PyCoupling;
    MSparse          Laplacian, EigenMatrix, Identity;
    vector<MatrixXf> FullEigenVectors, SortedEigenVectors;
    vector<VectorXf> FullEigenValues, SortedEigenValues;
    VectorXf         MeshGradient;

  EigenSolving(ndarray& Mesh,
               ndarray& PyMeshGradient,
               size_t   nMode,
               size_t   MaxIter,
               float    Tolerance,
               bool     debug=false): BaseLaplacian(Mesh){
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

                 this->debug             = debug;

                 float *adress           = (float*) PyMeshGradient.request().ptr;

                 Map<VectorXf> MeshGradient( adress, Nx * Ny );
                 this->MeshGradient = MeshGradient;
               }




   void Setlambda(float value){ this->lambda = value; this->k = 2.0 * PI / lambda; }

   float Getdx(){ return dx; }

   float Getdy(){ return dy; }

   float Getlambda(){ return lambda;}

   MSparse ComputeMatrix();

   tuple<MatrixXf, VectorXf> ComputeEigen(float alpha);

   void LoopOverITR(ndarray ITRList, float alpha, size_t order);

   tuple<ndarray, ndarray> GetSlice(uint slice);

   ndarray ComputingOverlap();

   Cndarray ComputingCoupling();

   Cndarray ComputingAdiabatic();

   ndarray ComputingIndices();

   void SortModesFields();

   void SortModesIndex();

   void ComputeLaplacian();

   void PringLaplacian();
};
