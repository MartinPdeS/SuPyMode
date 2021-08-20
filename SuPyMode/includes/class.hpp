class BaseLaplacian{
  public:
    int        TopSymmetry, BottomSymmetry, LeftSymmetry, RightSymmetry, Order;
    ndarray    Mesh;
    size_t     Nx, Ny, size;
    ScalarType dx, dy, D0xy, D1y, D2y, D1x, D2x;
    MSparse    Laplacian;
    bool       Debug;

    BaseLaplacian(ndarray&  Mesh, bool Debug){
      this->Nx                = Mesh.request().shape[0];
      this->Ny                = Mesh.request().shape[1];
      this->size              = Mesh.request().size;
      this->TopSymmetry      = 0;
      this->BottomSymmetry   = 0;
      this->LeftSymmetry     = 0;
      this->RightSymmetry    = 0;
      this->Debug            = Debug;
    }

    void SetLeftSymmetry(int value);
    int GetLeftSymmetry(){ return LeftSymmetry; }

    void SetRightSymmetry(int value);
    int GetRightSymmetry(){ return RightSymmetry; }

    void SetTopSymmetry(int value);
    int GetTopSymmetry(){ return TopSymmetry; }

    void SetBottomSymmetry(int value);
    int GetBottomSymmetry(){ return BottomSymmetry; }


    void SetLeftSymmetry3(int value);
    void SetRightSymmetry3(int value);
    void SetTopSymmetry3(int value);
    void SetBottomSymmetry3(int value);

    void SetLeftSymmetry5(int value);
    void SetRightSymmetry5(int value);
    void SetTopSymmetry5(int value);
    void SetBottomSymmetry5(int value);

    void Points3Laplacian();

    void Laplacian3Boundary();

    void Laplacian5Boundary();

    MSparse Points5Laplacian();

    void Setdx(ScalarType value){ dx = value; }

    void Setdy(ScalarType value){ dy = value; }
};


class EigenSolving : public BaseLaplacian{
  public:
    size_t             nMode, sMode, MaxIter, Nx, Ny, size, DegenerateFactor, ITRLength, Order, ExtrapolOrder;
    ScalarType         Tolerance, dx, dy, k, kInit, kDual, lambda, lambdaInit, MaxIndex, alpha;
    ScalarType        *MeshPtr, *ITRPtr;
    ndarray            Mesh, ITRList, PyOverlap, PyIndices;
    MSparse            Laplacian, EigenMatrix, Identity, M;
    vector<MatrixType> FullEigenVectors, SortedEigenVectors;
    vector<VectorType> FullEigenValues, SortedEigenValues;
    VectorType         MeshGradient;
    bool               Debug;
    BiCGSTAB<MSparse>  Solver;

  EigenSolving(ndarray&   Mesh,
               ndarray&   PyMeshGradient,
               size_t     nMode,
               size_t     sMode,
               size_t     MaxIter,
               ScalarType Tolerance,
               bool       Debug): BaseLaplacian(Mesh, Debug){
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
                 this->Debug             = Debug;
                 ScalarType *adress           = (ScalarType*) PyMeshGradient.request().ptr;

                 Map<VectorType> MeshGradient( adress, size );
                 this->MeshGradient = MeshGradient;

                 ComputeMaxIndex();


               }




   void Setlambda(ScalarType value){ this->lambda = value; this->k = 2.0 * PI / lambda;}


   void    LoopOverITR(ndarray ITRList, size_t order);

   void    LoopOverITR_(ndarray ITRList, size_t order);


   ScalarType              Getdx(){ return dx; }

   ScalarType              Getdy(){ return dy; }

   ScalarType              Getlambda(){ return lambda;}

   tuple<ndarray, ndarray> GetSlice(size_t slice);

   ndarray                 GetFields(size_t slice);

   ndarray                 GetIndices();

   ndarray                 GetBetas();


   tuple<MatrixType, VectorType> ComputeEigen(ScalarType alpha);

   ScalarType                    ComputeMaxIndex();

   MSparse                       ComputeMatrix();

   tuple<VectorType, ScalarType> ComputePreSolution(size_t& slice, size_t& mode);

   ndarray                       ComputingOverlap();

   ndarray                       ComputingCoupling();

   ndarray                       ComputingAdiabatic();

   vector<VectorType>            ComputeBetas();

   vector<size_t>                ComputecOverlaps(MatrixType Matrix0, MatrixType Matrix1, size_t idx);

   void                          ComputePSMMatrix();

   void                          ComputeLaplacian(size_t order);



   void SortModesFields();

   void SortModesIndex();

   void PSM(VectorType& X0, size_t& slice, size_t& mode);
};
