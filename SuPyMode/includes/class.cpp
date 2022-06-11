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

}


void
EigenSolving::PopulateModes(size_t Slice, MatrixType& EigenVectors, VectorType& EigenValues)
{
  for (SuperMode& mode : SuperModes)
  {
    mode.Fields.col(Slice)   << EigenVectors.col(mode.ModeNumber);
    mode.Fields.col(Slice).normalize();
    mode.Betas[Slice]        = sqrt( - EigenValues[mode.ModeNumber] ) / ITRList[Slice];
    mode.EigenValues[Slice]  = EigenValues[mode.ModeNumber];
    mode.Index[Slice]        = sqrt( abs( mode.EigenValues[Slice] ) ) / (ITRList[Slice] * kInit);
  }

}


void
EigenSolving::LoopOverITR(std::vector<double> ITRList, size_t order = 1){
  this->ITRList     = ITRList;

  ITRLength        = ITRList.size();


  kInit            = 2.0 * PI / lambda;

  MatrixType EigenVectors;
  VectorType EigenValues;

  PrepareSuperModes();


  std::vector<ScalarType> AllFirstEigenValues;
  AllFirstEigenValues.reserve(ITRLength);


  ScalarType alpha = -pow( k * ComputeMaxIndex(), 2 );



  for (size_t slice=0; slice<ITRLength; ++slice)
  {

    if (Debug)
    {
      size_t barWidth = 70;
      std::cout << "[";

      double progress = (double) slice/ITRLength;

      size_t pos = (size_t) (barWidth * progress);

      for (size_t i = 0; i < barWidth; ++i) {
          if (i < pos) std::cout << "=";
          else if (i == pos) std::cout << ">";
          else std::cout << " ";
      }
      std::cout << "] " << "ITR: " <<slice << "\n";
      std::cout.flush();

    }





    kDual = kInit * ITRList[slice];

    tie(EigenVectors, EigenValues) = ComputeEigen(alpha);

    PopulateModes(slice, EigenVectors, EigenValues);

    AllFirstEigenValues.push_back(EigenValues[0]);

    size_t next = slice+1, mode=0;

    alpha = ExtrapolateNext(order, AllFirstEigenValues, ITRList, next);
  }

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
              std::cout<<"Bad mode correspondence: "<< BestOverlap <<"  At ITR: "<< ITRList[Slice] <<". You should consider makes more ITR steps"<<std::endl;
    }

  return Indices;

}




void
EigenSolving::SortModes(std::string Type)
{
  std::cout<<"Sorting SuperModes\n";
  SortedSuperModes = std::vector<SuperMode>(sMode);

  for (SuperMode &mode : SortedSuperModes)
      mode.Init(ITRLength, Nx, Ny, LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry, sMode);

  if (Type == "Field") SortModesFields();
  else if (Type == "Index") SortModesIndex();
  else if (Type == "None") SortModesNone();
}










void
EigenSolving::SortSliceIndex(size_t Slice)
{
  vector<ScalarType> Betas;
  Betas.reserve(nMode);


  size_t iter=0;
  for (size_t mode=0; mode<sMode; ++mode)
  {
      Betas.push_back(SuperModes[mode].Betas[Slice]);

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
EigenSolving::SortModesFields()
{
  for (size_t mode=0; mode<sMode; ++mode)
      SortedSuperModes[mode] = SuperModes[mode];


  for (size_t slice=0; slice<ITRLength-1; ++slice)
      SortSliceFields(slice);
}


void
EigenSolving::SortSliceFields(size_t Slice)
{
  for (size_t previous=0; previous<sMode; ++previous)
  {
      SuperMode &Mode0 = SortedSuperModes[previous];
      std::vector<ScalarType> Overlaps(nMode, 0);

      for (size_t after=0; after<nMode; ++after)
          {
            SuperMode &Mode1 = SuperModes[after];
            Overlaps[after] = abs( Mode0.Fields.col(Slice).transpose() * Mode1.Fields.col(Slice+1)  );
          }

    int bestFit = std::max_element(Overlaps.begin(), Overlaps.end()) - Overlaps.begin();

    Mode0.Fields.col(Slice+1) = SuperModes[bestFit].Fields.col(Slice+1);
    Mode0.Betas[Slice+1]      = SuperModes[bestFit].Betas[Slice+1];
    Mode0.Index[Slice+1]      = SuperModes[bestFit].Index[Slice+1];
  }
}



void
EigenSolving::SortModesNone()
{
  SortedSuperModes = std::vector<SuperMode>(sMode);

  for (size_t mode=0; mode<sMode; ++mode)
      SortedSuperModes[mode] = SuperModes[mode];

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













// -
