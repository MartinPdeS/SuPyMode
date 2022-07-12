#pragma once


#include "Definitions.cpp"

struct SuperMode
{
  VectorType Betas, EigenValues, Index;
  size_t ITRLength, Nx, Ny, ModeNumber, sMode;
  int LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry;
  MatrixType Fields, Coupling, Adiabatic;

  SuperMode(size_t ModeNumber){this->ModeNumber = ModeNumber;}
  SuperMode(){}


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
      this->Fields.col(Slice) = Other.Fields.col(Slice);
      this->Betas[Slice]      = Other.Betas[Slice];
      this->Index[Slice]      = Other.Index[Slice];
      this->Adiabatic         = Other.Adiabatic;
      this->Coupling          = Other.Coupling;
  }


  ScalarType ComputeOverlap(SuperMode& Other, size_t Slice)
  {
    return 1;//abs( this->Fields.col(Slice) * Other.Fields.col(Slice) );
  }


  ScalarType ComputeCoupling(SuperMode& Other, size_t Slice, VectorType &MeshGradient, ScalarType &kInit)
  {
    ComplexScalarType C;
    if (this->ModeNumber == Other.ModeNumber)
    {
      C = 0.0;
    }

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
    ComplexScalarType C;
    ScalarType A;

    ScalarType Beta0 = this->Betas[Slice], Beta1 = Other.Betas[Slice];

    if (this->ModeNumber == Other.ModeNumber)
    {
      C = 0.0;
      A = 0.0;
    }

    else
    {
      VectorType overlap = this->Fields.col(Slice).cwiseProduct( Other.Fields.col(Slice) );

      ScalarType Beta0 = this->Betas[Slice], Beta1 = Other.Betas[Slice];

      C  = - (ScalarType) 0.5 * J * kInit*kInit / sqrt(Beta0 *  Beta1) * abs( 1.0f / (Beta0 - Beta1) );

      ScalarType I       = Trapz(overlap.cwiseProduct( MeshGradient ), 1.0, Nx, Ny);

      C      *=  I;

      C = abs(C);

      A = abs( (Beta0-Beta1) / C );
    }

    this->Coupling(Other.ModeNumber, Slice) = abs(C);
    Other.Coupling(this->ModeNumber, Slice) = abs(C);
    this->Adiabatic(Other.ModeNumber, Slice) = A;
    Other.Adiabatic(this->ModeNumber, Slice) = A;
  }


  ndarray GetFields(){ return Eigen2ndarray_( this->Fields, { ITRLength, Nx, Ny} ); }
  ndarray GetIndex(){ return Eigen2ndarray_( this->Index, { ITRLength} ); }
  ndarray GetBetas(){ return Eigen2ndarray_( this->Betas, { ITRLength} ); }
  ndarray GetAdiabatic(){ return Eigen2ndarray_( this->Adiabatic, { ITRLength, sMode} ); }

  ndarray GetCoupling(){ return Eigen2ndarray_( this->Coupling, { ITRLength, sMode} ); }
  ndarray GetAdiabaticSpecific(SuperMode& Mode){ return Eigen2ndarray_( this->Adiabatic.row(Mode.ModeNumber), { ITRLength} ); }
  ndarray GetCouplingSpecific(SuperMode& Mode){ return Eigen2ndarray_( this->Coupling.row(Mode.ModeNumber), { ITRLength} ); }
};
