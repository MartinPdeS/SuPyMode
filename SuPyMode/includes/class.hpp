std::complex<ScalarType> J(0.0, 1.0);


struct SuperMode
{
  MatrixType Fields;
  VectorType Betas, EigenValues;
  size_t ITRLength, Nx, Ny, ModeNumber;

  SuperMode(size_t ModeNumber){this->ModeNumber = ModeNumber;}
  SuperMode(){}

  ndarray GetField(size_t slice)
  {
    MatrixType * Vectors = new MatrixType;

    (*Vectors) = this->Fields.col(slice);

    return Eigen2ndarray( Vectors, { Nx, Ny}  );
  }

  void Init(size_t ITRLength, size_t Nx, size_t Ny)
  {
    this->Nx = Nx;
    this->Ny = Ny;
    this->ITRLength   = ITRLength;
    this->Fields      = MatrixType(Nx * Ny, ITRLength);
    this->Betas       = VectorType(ITRLength);
    this->EigenValues = VectorType(ITRLength);
  }

  void CopyOtherSlice(SuperMode& Other, size_t Slice)
  {
      Fields.col(Slice) = Other.Fields.col(Slice);
      Betas[Slice]  = Other.Betas[Slice];
  }


  ScalarType ComputeOverlap(SuperMode& Other, size_t Slice)
  {
    return abs( this->Fields.col(Slice).transpose() * Other.Fields.col(Slice) );
  }


  ScalarType ComputeCoupling(SuperMode& Other, size_t Slice, VectorType &MeshGradient, ScalarType &kInit)
  {
    VectorType overlap = this->Fields.col(Slice).cwiseProduct( Other.Fields.col(Slice) );

    ScalarType Beta0 = this->Betas[Slice],
               Beta1 = Other.Betas[Slice];


    ComplexScalarType C  = - (ScalarType) 0.5 * J * kInit*kInit / sqrt(Beta0 *  Beta1) * abs( 1.0f / (Beta0 - Beta1) );

    ScalarType I       = Trapz(overlap.cwiseProduct( MeshGradient ), 1.0, Nx, Ny);

    C      *=  I;

    return abs(C);

  }



  ndarray GetFields()
  {
    MatrixType * Vectors = new MatrixType;

    (*Vectors) = this->Fields;

    return Eigen2ndarray( Vectors, { ITRLength, Nx, Ny} );
  }

};




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
    MatrixType         Mode0, Mode1, Mode2, Mode3;
    SuperMode          SuperMode0, SuperMode1, SuperMode2;
    std::vector<SuperMode>        SuperModes;

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
               : BaseLaplacian(Mesh, Debug)
                {
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

                 for (int i=0; i<nMode; ++i)
                    SuperModes.push_back(SuperMode(i));

                 ComputeMaxIndex();


               }



   void PopulateModes(size_t Slice, MatrixType& EigenVectors, VectorType& EigenValues);
   ndarray GetMode(size_t Mode);
   void PrepareSuperModes();
   void SwapMode(SuperMode &Mode0, SuperMode &Mode1);
   void SortSliceIndex(size_t Slice);

   void    LoopOverITR(ndarray ITRList, size_t order);

   tuple<ndarray, ndarray> GetSlice(size_t slice);

   ndarray                 GetFullEigen();

   ndarray                 GetFields(size_t slice);

   ndarray                 GetIndices();

   ndarray                 GetBetas();

   tuple<MatrixType, VectorType> ComputeEigen(ScalarType alpha);

   ScalarType                    ComputeMaxIndex();

   MSparse                       ComputeMatrix();

   ndarray                       ComputingOverlap();

   ndarray                       ComputingCoupling();

   ndarray                       ComputingAdiabatic();

   vector<VectorType>            ComputeBetas();

   vector<size_t>                ComputecOverlaps(MatrixType Matrix0, MatrixType Matrix1, size_t idx);

   void                          ComputeLaplacian(size_t order);



   void SortModesFields();

   void SortModesIndex();

   void SortModes(std::string Type);

};
