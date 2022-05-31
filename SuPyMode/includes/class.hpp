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
      this->Debug             = Debug;
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
};


class EigenSolving : public BaseLaplacian{
  public:
    size_t             nMode, sMode, MaxIter, Nx, Ny, size, DegenerateFactor, ITRLength, Order, ExtrapolOrder;
    ScalarType         Tolerance, k, kInit, kDual, lambda, lambdaInit, MaxIndex, alpha;
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
               ScalarType dx,
               ScalarType dy,
               ScalarType Wavelength,
               bool       Debug)
               : BaseLaplacian(Mesh, Debug){
                 this->dx                = dx;
                 this->dy                = dy;
                 this->nMode             = nMode;
                 this->sMode             = sMode;
                 this->MaxIter           = MaxIter;
                 this->Tolerance         = Tolerance;
                 this->Mesh              = Mesh;
                 this->Nx                = Mesh.request().shape[0];
                 this->Ny                = Mesh.request().shape[1];
                 this->size              = Mesh.request().size;
                 this->MeshPtr           = (ScalarType*) Mesh.request().ptr;
                 this->Debug             = Debug;
                 this->lambda            = Wavelength;
                 this->k                 = 2.0 * PI / Wavelength;
                 this->kInit             = this->k;
                 ScalarType *adress      = (ScalarType*) PyMeshGradient.request().ptr;

                 Map<VectorType> MeshGradient( adress, size );
                 this->MeshGradient = MeshGradient;

                 ComputeMaxIndex();


               }



   void    LoopOverITR(ndarray ITRList, size_t order);

   void    LoopOverITR_(ndarray ITRList, size_t order);

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

   void SortModes(std::string Type);

   void PSM(VectorType& X0, size_t& slice, size_t& mode);
};
