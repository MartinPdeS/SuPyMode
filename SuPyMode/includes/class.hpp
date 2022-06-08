std::complex<ScalarType> J(0.0, 1.0);


struct SuperMode
{
  MatrixType Fields;
  VectorType Betas, EigenValues, Index;
  size_t ITRLength, Nx, Ny, ModeNumber, sMode;
  int LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry;
  MatrixType Coupling, Adiabatic;

  SuperMode(size_t ModeNumber){this->ModeNumber = ModeNumber;}
  SuperMode(){}

  ndarray GetField(size_t slice)
  {
    MatrixType * Vectors = new MatrixType;

    (*Vectors) = this->Fields.col(slice);

    return Eigen2ndarray( Vectors, { Nx, Ny}  );
  }

  void Init(size_t &ITRLength,
            size_t &Nx,
            size_t &Ny,
            int    &LeftSymmetry,
            int    &RightSymmetry,
            int    &TopSymmetry,
            int    &BottomSymmetry,
            int     sMode)
  {
    this->Nx             = Nx;
    this->Ny             = Ny;
    this->sMode          = sMode;
    this->ITRLength      = ITRLength;
    this->Fields         = MatrixType(Nx * Ny, ITRLength);
    this->Betas          = VectorType(ITRLength);
    this->EigenValues    = VectorType(ITRLength);
    this->Index          = VectorType(ITRLength);
    this->Adiabatic      = MatrixType(sMode, ITRLength);
    this->Coupling       = MatrixType(sMode, ITRLength);
    this->BottomSymmetry = BottomSymmetry;
    this->TopSymmetry    = TopSymmetry;
    this->RightSymmetry  = RightSymmetry;
    this->LeftSymmetry   = LeftSymmetry;
  }

  void CopyOtherSlice(SuperMode& Other, size_t Slice)
  {
      Fields.col(Slice) = Other.Fields.col(Slice);
      Betas[Slice]      = Other.Betas[Slice];
      Index[Slice]      = Other.Index[Slice];
      Adiabatic         = Other.Adiabatic;
      Coupling          = Other.Coupling;
  }


  ScalarType ComputeOverlap(SuperMode& Other, size_t Slice)
  {
    return abs( this->Fields.col(Slice).transpose() * Other.Fields.col(Slice) );
  }


  ScalarType ComputeCoupling(SuperMode& Other, size_t Slice, VectorType &MeshGradient, ScalarType &kInit)
  {
    ComplexScalarType C;
    if (this->ModeNumber == Other.ModeNumber){C = 0.0;}

    else
    {
      VectorType overlap = this->Fields.col(Slice).cwiseProduct( Other.Fields.col(Slice) );

      ScalarType Beta0 = this->Betas[Slice], Beta1 = Other.Betas[Slice];

      C  = - (ScalarType) 0.5 * J * kInit*kInit / sqrt(Beta0 *  Beta1) * abs( 1.0f / (Beta0 - Beta1) );

      ScalarType I       = Trapz(overlap.cwiseProduct( MeshGradient ), 1.0, Nx, Ny);

      C      *=  I;

      C = abs(C);
    }

    this->Coupling(Other.ModeNumber, Slice) = abs(C);
    Other.Coupling(this->ModeNumber, Slice) = abs(C);

    return abs(C);
  }


  ScalarType ComputeAdiabatic(SuperMode& Other, size_t Slice, VectorType &MeshGradient, ScalarType &kInit)
  {
    ScalarType A;

    ScalarType Beta0 = this->Betas[Slice], Beta1 = Other.Betas[Slice];

    if (this->ModeNumber == Other.ModeNumber) { A = 0.0; }
    else { A = abs( (Beta0-Beta1) / ComputeCoupling(Other, Slice, MeshGradient, kInit) ); }

    this->Adiabatic(Other.ModeNumber, Slice) = A;
    Other.Adiabatic(this->ModeNumber, Slice) = A;

    return A;
  }


  void PopulateCouplingAdiabatic(SuperMode& Other, size_t Slice, VectorType &MeshGradient, ScalarType &kInit)
  {
    ComputeCoupling(Other, Slice, MeshGradient, kInit);
    ComputeAdiabatic(Other, Slice, MeshGradient, kInit);
  }






  ndarray GetFields(){ return Eigen2ndarray_( this->Fields, { ITRLength, Nx, Ny} ); }
  ndarray GetIndex(){ return Eigen2ndarray_( this->Index, { ITRLength} ); }
  ndarray GetBetas(){ return Eigen2ndarray_( this->Betas, { ITRLength} ); }
  ndarray GetAdiabatic(){ return Eigen2ndarray_( this->Adiabatic, { ITRLength, sMode} ); }
  ndarray GetCoupling(){ return Eigen2ndarray_( this->Coupling, { ITRLength, sMode} ); }

};






class BaseLaplacian{
  public:
    int        TopSymmetry, BottomSymmetry, LeftSymmetry, RightSymmetry, Order;
    ndarray    Mesh;
    size_t     Nx, Ny, size;
    ScalarType dx, dy, D0xy, D1y, D2y, D1x, D2x;
    MSparse    Laplacian;
    bool       Debug;

    BaseLaplacian(ndarray&  Mesh, ScalarType dx, ScalarType dy){
      this->Nx                = Mesh.request().shape[0];
      this->Ny                = Mesh.request().shape[1];
      this->size              = Mesh.request().size;
      this->dx                = dx;
      this->dy                = dy;
      this->Mesh              = Mesh;
    }

    void SetLeftSymmetry();
    void SetRightSymmetry();
    void SetTopSymmetry();
    void SetBottomSymmetry();

    void SetLeftSymmetry3();
    void SetRightSymmetry3();
    void SetTopSymmetry3();
    void SetBottomSymmetry3();

    void SetLeftSymmetry5();
    void SetRightSymmetry5();
    void SetTopSymmetry5();
    void SetBottomSymmetry5();

    void Points3Laplacian();

    void Laplacian3Boundary();

    void Laplacian5Boundary();

    MSparse Points5Laplacian();
};


class EigenSolving : public BaseLaplacian
{
  public:
    size_t             nMode, sMode, MaxIter, DegenerateFactor, ITRLength, Order, ExtrapolOrder;
    ScalarType         Tolerance, k, kInit, kDual, lambda, MaxIndex;
    ScalarType        *MeshPtr, *ITRPtr;
    ndarray            ITRList;
    MSparse            EigenMatrix, Identity, M;
    VectorType         MeshGradient;
    BiCGSTAB<MSparse>  Solver;
    std::vector<SuperMode> SuperModes, SortedSuperModes;

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
               : BaseLaplacian(Mesh, dx, dy)
                {
                 this->Debug             = Debug;
                 this->nMode             = nMode;
                 this->sMode             = sMode;
                 this->MaxIter           = MaxIter;
                 this->Tolerance         = Tolerance;

                 this->MeshPtr           = (ScalarType*) Mesh.request().ptr;
                 this->lambda            = Wavelength;
                 this->k                 = 2.0 * PI / Wavelength;
                 this->kInit             = this->k;
                 ScalarType *adress      = (ScalarType*) PyMeshGradient.request().ptr;

                 Map<VectorType> MeshGradient( adress, size );
                 this->MeshGradient = MeshGradient;

                 for (int i=0; i<nMode; ++i)
                    SuperModes.push_back(SuperMode(i));

                 for (int i=0; i<sMode; ++i)
                    SortedSuperModes.push_back(SuperMode(i));

                 ComputeMaxIndex();


               }



   void      PopulateModes(size_t Slice, MatrixType& EigenVectors, VectorType& EigenValues);
   SuperMode GetMode(size_t Mode){ return SortedSuperModes[Mode]; }

   void      PrepareSuperModes();
   void      SwapMode(SuperMode &Mode0, SuperMode &Mode1);
   void      SortSliceIndex(size_t Slice);

   void      LoopOverITR(ndarray ITRList, size_t order);

   tuple<ndarray, ndarray> GetSlice(size_t slice);
   ndarray                 GetIndices();
   ndarray                 GetBetas();
   ndarray                 GetFields();

   tuple<MatrixType, VectorType> ComputeEigen(ScalarType alpha);

   ScalarType                    ComputeMaxIndex();

   MSparse                       ComputeMatrix();

   void                          ComputeCoupling();
   void                          ComputeAdiabatic();
   vector<size_t>                ComputecOverlaps(size_t idx);
   void                          ComputeLaplacian(size_t order);
   void                          ComputeCouplingAdiabatic();

   void SortModesFields();
   void SortModesIndex();
   void SortModes(std::string Type);

};
