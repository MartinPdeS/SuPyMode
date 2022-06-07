#include "class.hpp"
#include "Extrapolate.hpp"
#include "Laplacian.cpp"




MSparse
EigenSolving::ComputeMatrix(){

    EigenMatrix = BaseLaplacian::Laplacian;

    size_t iter = 0;

    for(size_t i=0; i<Nx; ++i)
       for(size_t j=0; j<Ny; ++j){
           Identity.coeffRef(iter,iter) = + pow( MeshPtr[iter] * kDual, 2);
           ++iter;
         }

    EigenMatrix += Identity;

    Identity.setIdentity();

    return -1.0*EigenMatrix;
    }


tuple<MatrixType, VectorType>
EigenSolving::ComputeEigen(ScalarType alpha){

    MSparse EigenMatrix = ComputeMatrix();

    SparseGenRealShiftSolve<ScalarType> op(EigenMatrix);

    GenEigsRealShiftSolver<SparseGenRealShiftSolve<ScalarType>> eigs(op, nMode, 2*nMode, alpha);

    eigs.init();

    int nconv = eigs.compute(SortRule::LargestMagn, MaxIter, Tolerance);

    MatrixType Vectors = eigs.eigenvectors().real();

    VectorType Values = eigs.eigenvalues().real();

    return std::make_tuple(Vectors, Values);
    }


void
EigenSolving::ComputeLaplacian(size_t Order){
  Identity = MSparse(size,size); Identity.setIdentity();

  this->Order = Order;

  switch(Order){

  case 2:
      Points3Laplacian();

  case 4:
      Points5Laplacian();

  SetLeftSymmetry();
  SetRightSymmetry();
  SetTopSymmetry();
  SetBottomSymmetry();
  }

}
























void
EigenSolving::PrepareSuperModes()
{
  for (SuperMode& mode : SuperModes)
      mode.Init(ITRLength, Nx, Ny, LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry, sMode);

  for (SuperMode& mode : SortedSuperModes)
      mode.Init(ITRLength, Nx, Ny, LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry, sMode);
}


void
EigenSolving::PopulateModes(size_t Slice, MatrixType& EigenVectors, VectorType& EigenValues)
{
  for (SuperMode& mode : SuperModes)
  {
    mode.Fields.col(Slice)   << EigenVectors.col(mode.ModeNumber);
    mode.Betas[Slice]        = sqrt( - EigenValues[mode.ModeNumber] ) / ITRPtr[Slice];
    mode.EigenValues[Slice]  = EigenValues[mode.ModeNumber];
    mode.Index[Slice]        = sqrt( abs( mode.EigenValues[Slice] ) ) / (ITRPtr[Slice] * kInit);
  }

}


void
EigenSolving::LoopOverITR(ndarray ITRList, size_t order = 1){


  ITRList          = ITRList;
  ITRLength        = ITRList.request().size;
  ITRPtr           = (ScalarType*) ITRList.request().ptr;

  kInit            = 2.0 * PI / lambda;

  MatrixType EigenVectors;
  VectorType EigenValues;

  PrepareSuperModes();


  std::vector<ScalarType> AllFirstEigenValues;
  AllFirstEigenValues.reserve(ITRLength);


  ScalarType alpha = -pow( k * ComputeMaxIndex(), 2 );

  size_t barWidth = 70;
  std::cout << "[";
  for (size_t slice=0; slice<ITRLength; ++slice)
  {
    double progress = (double) slice/ITRLength;

    size_t pos = (size_t) (barWidth * progress);

    for (size_t i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << "ITR: " <<slice << "\n";
    std::cout.flush();



    kDual = kInit * ITRPtr[slice] ;

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);

    PopulateModes(slice, EigenVectors, EigenValues);

    AllFirstEigenValues.push_back(EigenValues[0]);

    size_t next = slice+1, mode=0;

    alpha = ExtrapolateNext(order, AllFirstEigenValues, ITRList, next);
  }
  SortModesIndex();
}




































vector<size_t>
EigenSolving::ComputecOverlaps(size_t Slice){

  ScalarType BestOverlap, Overlap;
  vector<size_t> Indices(sMode);

  for (size_t i=0; i<sMode; ++i)
    {
      BestOverlap = 0;
      for (size_t j=0; j<nMode; ++j)
          {
              SuperMode Mode0 = SuperModes[i], Mode1 = SuperModes[j];

              Overlap = Mode0.ComputeOverlap(Mode1, Slice);

              if (Overlap > BestOverlap) {Indices[i] = j; BestOverlap = Overlap;}
          }
          if (BestOverlap<0.8)
              std::cout<<"Bad mode correspondence: "<< BestOverlap <<"  At ITR: "<< ITRPtr[Slice] <<". You should consider makes more ITR steps"<<std::endl;
    }

  return Indices;

}


void
EigenSolving::SortModes(std::string Type)
{
  std::cout<<"Sorting SuperModes\n";
  if (Type == "Field") SortModesFields();
  else if (Type == "Index") SortModesIndex();
}


void
EigenSolving::SortModesFields(){

}


void
EigenSolving::SwapMode(SuperMode &Mode0, SuperMode &Mode1)
{
  MatrixType TempField0 = Mode0.Fields,
             TempField1 = Mode1.Fields;

  VectorType TempBeta0 = Mode0.Betas,
             TempBeta1 = Mode1.Betas;

  Mode0.Fields = TempField1;
  Mode1.Fields = TempField0;

  Mode0.Betas = TempBeta1;
  Mode1.Betas = TempBeta0;
}


void
EigenSolving::SortSliceIndex(size_t Slice)
{
  vector<ScalarType> Betas;
  Betas.reserve(nMode);
  SortedSuperModes = std::vector<SuperMode>(sMode);

  size_t iter=0;
  for (size_t mode=0; mode<sMode; ++mode)
  {
      Betas.push_back(SuperModes[mode].Betas[Slice]);
      SortedSuperModes[mode] = SuperModes[mode];
      ++iter;
  }

  vector<size_t> sorted = sort_indexes( Betas );

  for (size_t mode=0; mode<sMode; ++mode)
  {
      auto order = sorted[mode];

      SortedSuperModes[mode].CopyOtherSlice(SuperModes[order], Slice);
  }

}


void
EigenSolving::SortModesIndex()
{
  for (size_t l=0; l<ITRLength; ++l)
      SortSliceIndex(l);

}

void
EigenSolving::ComputeCoupling()
{

  std::cout<<"Computing coupling\n";

  for (size_t slice=0; slice<ITRLength; ++slice)
      for (SuperMode &mode0 : SortedSuperModes)
          for (size_t m=0; m<sMode;++m)
              {
                SuperMode &mode1 = SortedSuperModes[m];

                mode1.ComputeCoupling(mode0, slice, MeshGradient, kInit);

              }

}


void
EigenSolving::ComputeAdiabatic(){

  std::cout<<"Computing adiabatic\n";

  for (size_t slice=0; slice<ITRLength; ++slice)
      for (SuperMode &mode0 : SortedSuperModes)
          for (size_t m=0; m<sMode;++m)
              {
                SuperMode &mode1 = SortedSuperModes[m];

                mode1.Adiabatic(mode0.ModeNumber, slice) = mode1.ComputeAdiabatic(mode0, slice, MeshGradient, kInit);
              }

}




void
EigenSolving::ComputeCouplingAdiabatic(){

  std::cout<<"Computing coupling/adiabatic\n";

  for (size_t slice=0; slice<ITRLength; ++slice)
      for (SuperMode &mode0 : SortedSuperModes)
          for (size_t m=0; m<sMode;++m)
                mode0.PopulateCouplingAdiabatic(SortedSuperModes[m], slice, MeshGradient, kInit);

}

ScalarType
EigenSolving::ComputeMaxIndex(){
  MaxIndex = 0.0;
  for (size_t i=0; i<size; ++i)
     if (MeshPtr[i] > MaxIndex)
         MaxIndex = MeshPtr[i];

 return MaxIndex;
}


tuple<ndarray, ndarray>
EigenSolving::GetSlice(size_t Slice){
  MatrixType OutputFields(size, sMode);
  VectorType OutputBetas(sMode);

  for (size_t mode=0; mode<sMode; ++mode)
  {
    OutputFields.col(mode) = SuperModes[mode].Fields.col(Slice);
    OutputBetas[mode]      = SuperModes[mode].Betas[Slice];
  }


  ndarray FieldsPython = Eigen2ndarray_( OutputFields, { sMode, Nx, Ny } ) ;

  ndarray BetasPython = Eigen2ndarray_( OutputBetas, { sMode } ) ;

  return std::make_tuple( FieldsPython, BetasPython );
}





ndarray
EigenSolving::GetFields(){
  MatrixType Output(ITRLength*size, sMode);

  for (size_t mode=0; mode<sMode; ++mode)
    Output.col(mode)  = SortedSuperModes[mode].Fields;

  return Eigen2ndarray_( Output, { sMode, ITRLength, Nx, Ny } ) ;
}


ndarray
EigenSolving::GetIndices(){
  MatrixType Output(ITRLength, sMode);

  for (size_t mode=0; mode<sMode; ++mode)
    Output.col(mode)  = SortedSuperModes[mode].Index;

  return Eigen2ndarray_( Output, { ITRLength, sMode } ) ;
}


ndarray
EigenSolving::GetBetas(){
  MatrixType Output(ITRLength, sMode);

  for (size_t mode=0; mode<sMode; ++mode)
    Output.col(mode)  = SortedSuperModes[mode].Betas;

  return Eigen2ndarray_( Output, { ITRLength, sMode } ) ;
}












// -
