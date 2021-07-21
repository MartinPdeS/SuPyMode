class BaseLaplacian{
  public:
    int      TopSymmetry, BottomSymmetry, LeftSymmetry, RightSymmetry;
    ndarray  Mesh;
    size_t   Nx, Ny, size;
    ScalarType    dx, dy;
    MSparse  Laplacian;

    BaseLaplacian(ndarray&  Mesh){
      this->Nx                = Mesh.request().shape[0];
      this->Ny                = Mesh.request().shape[1];
      this->size              = Mesh.request().size;
      this->TopSymmetry      = 0;
      this->BottomSymmetry   = 0;
      this->LeftSymmetry     = 0;
      this->RightSymmetry    = 0;
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

    MSparse Points5Laplacian();

    void Setdx(ScalarType value){ dx = value; }

    void Setdy(ScalarType value){ dy = value; }
};


class EigenSolving : public BaseLaplacian{
  public:
    size_t             nMode, sMode, MaxIter, Nx, Ny, size, DegenerateFactor, ITRLength;
    ScalarType         Tolerance, dx, dy, k, kInit, kDual, lambda, lambdaInit, MaxIndex;
    ScalarType        *MeshPtr, *ITRPtr;
    ndarray            Mesh, ITRList, PyOverlap, PyIndices, PyFullEigenVectors, PyFullEigenValues, PyAdiabatic;
    Cndarray           PyCoupling;
    MSparse            Laplacian, EigenMatrix, Identity;
    vector<MatrixType> FullEigenVectors, SortedEigenVectors;
    vector<VectorType> FullEigenValues, SortedEigenValues;
    VectorType         MeshGradient;

  EigenSolving(ndarray& Mesh,
               ndarray& PyMeshGradient,
               size_t   nMode,
               size_t   sMode,
               size_t   MaxIter,
               ScalarType    Tolerance): BaseLaplacian(Mesh){
                 this->nMode             = nMode;
                 this->sMode             = sMode;
                 this->MaxIter           = MaxIter;
                 this->Tolerance         = Tolerance;
                 this->Mesh              = Mesh;
                 this->Nx                = Mesh.request().shape[0];
                 this->Ny                = Mesh.request().shape[1];
                 this->size              = Mesh.request().size;
                 this->dx                = 100./this->Nx;
                 this->dy                = 100./this->Ny;
                 this->MeshPtr           = (ScalarType*) Mesh.request().ptr;
                 this->DegenerateFactor  = 1.0;

                 ScalarType *adress           = (ScalarType*) PyMeshGradient.request().ptr;

                 Map<VectorType> MeshGradient( adress, size );
                 this->MeshGradient = MeshGradient;

                 ComputeMaxIndex();


               }


   ScalarType ComputeMaxIndex();

   void ComputeDegenerateFactor();

   void Setlambda(ScalarType value){ this->lambda = value; this->k = 2.0 * PI / lambda; }

   ScalarType Getdx(){ return dx; }

   ScalarType Getdy(){ return dy; }

   ScalarType Getlambda(){ return lambda;}

   MSparse ComputeMatrix();

   tuple<MatrixType, VectorType> ComputeEigen(ScalarType alpha);

   void LoopOverITR(ndarray ITRList, size_t order);

   tuple<ndarray, ndarray> GetSlice(size_t slice);

   ndarray ComputingOverlap();

   Cndarray ComputingCoupling();

   ndarray GetFields();

   ndarray GetIndices();

   ndarray GetBetas();

   Cndarray ComputingAdiabatic();

   void SortModesFields();

   void SortModesIndex();

   void ComputeLaplacian();

   void PringLaplacian();

   vector<VectorType> ComputeBetas();
};
