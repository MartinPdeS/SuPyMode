void
BaseLaplacian::Laplacian3Boundary(){

  for (size_t j=1; j<size-1; ++j){
      if (j%Ny == 0)
          Laplacian.coeffRef(j-1,j) = 0.0;

      if (j%Ny == Ny-1)
          Laplacian.coeffRef(j+1,j) = 0.0;
  }
}


void
BaseLaplacian::Laplacian5Boundary(){

  for (size_t j=1; j<size-1; ++j){
      if (j%Ny == Ny-1){
          Laplacian.coeffRef(j,j+1)   = 0.0;
          Laplacian.coeffRef(j,j+2)   = 0.0;
          Laplacian.coeffRef(j-1,j+1) = 0.0;
          ++j;
          Laplacian.coeffRef(j,j-1)   = 0.0;
          Laplacian.coeffRef(j,j-2)   = 0.0;
          Laplacian.coeffRef(j+1,j-1) = 0.0;
          }
  }
}

void
BaseLaplacian::Points3Laplacian(){
    Order = 2;
    D0xy  = SD2A2[1]*( 1./pow(dx,2) + 1./pow(dy,2) );
    D1y   = SD2A2[0]/pow(dy,2);
    D1x   = SD2A2[0]/pow(dx,2);

    MSparse Identity(size,size); Identity.setIdentity();
    Laplacian = MSparse(size,size);

    vector<T> Tri(size);
    for (size_t i=0; i<size-1; ++i){
        Tri.push_back(T(i+1, i, D1y));
        Tri.push_back(T(i, i+1, D1y));
        }

    for (size_t i=0; i<size-Ny; ++i){
        Tri.push_back(T(i+Ny, i, D1x));
        Tri.push_back(T(i, i+Ny, D1x));
        }

    Laplacian.setFromTriplets(Tri.begin(), Tri.end());

    Laplacian += Identity * D0xy;

    Laplacian3Boundary();

}



MSparse
BaseLaplacian::Points5Laplacian(){

  Order = 4;

  D0xy  = SD2A4[2]*( 1./pow(dx,2) + 1./pow(dy,2) );

  D1y  = SD2A4[1]/pow(dy,2);
  D2y  = SD2A4[0]/pow(dy,2);

  D1x  = SD2A4[1]/pow(dx,2);
  D2x  = SD2A4[0]/pow(dx,2);



  MSparse Identity(size,size); Identity.setIdentity();
  Laplacian = MSparse(size,size);

  vector<T> Tri(size);
  for (size_t i=0; i<size-1; ++i){
      Tri.push_back(T(i+1, i, D1y));
      Tri.push_back(T(i, i+1, D1y));
    }

  for (size_t i=0; i<size-2; ++i){
      Tri.push_back(T(i+2, i, D2y));
      Tri.push_back(T(i, i+2, D2y));
      }

  for (size_t i=0; i<size-Ny; ++i){
      Tri.push_back(T(i+Ny, i, D1x));
      Tri.push_back(T(i, i+Ny, D1x));
      }

  for (size_t i=0; i<size-2*Ny; ++i){
      Tri.push_back(T(i+2*Ny, i, D2x));
      Tri.push_back(T(i, i+2*Ny, D2x));
      }

  Laplacian.setFromTriplets(Tri.begin(), Tri.end());

  Laplacian += Identity * D0xy;

  Laplacian5Boundary();

  return Laplacian;
}






void BaseLaplacian::SetLeftSymmetry()
{
  if (LeftSymmetry == 0) return;

  std::cout<<"order: "<<Order<<"\n";
  switch (Order)
  {
    case 2: SetLeftSymmetry3();
    case 4: SetLeftSymmetry5();
  }
}

void BaseLaplacian::SetTopSymmetry()
{
  if (TopSymmetry == 0) return;
  switch (Order)
  {
    case 2: SetTopSymmetry3();
    case 4: SetTopSymmetry5();
  }
}

void BaseLaplacian::SetBottomSymmetry()
{
  if (BottomSymmetry == 0) return;
  switch (Order)
  {
    case 2: SetBottomSymmetry3();
    case 4: SetBottomSymmetry5();
  }
}

void BaseLaplacian::SetRightSymmetry()
{
  if (RightSymmetry == 0) return;
  switch (Order)
  {
    case 2: SetRightSymmetry3();
    case 4: SetRightSymmetry5();
  }
}




void BaseLaplacian::SetBottomSymmetry3()
{
  if (BottomSymmetry == 1)
      for (size_t j=0; j<size; ++j)
          if (j%Ny == 0)
              Laplacian.coeffRef(j,j+1)         += D1y;

  if (BottomSymmetry == -1)
      for (size_t j=0; j<size; ++j)
          if (j%Ny == 0)
              Laplacian.coeffRef(j,j+1)          = 0.0;
}


void BaseLaplacian::SetTopSymmetry3()
{
  if (TopSymmetry == 1)
      for (size_t j=1; j<size; ++j)
          if (j%Ny == Ny-1)
              Laplacian.coeffRef(j,j-1)          += D1y;


  if (TopSymmetry == -1)
      for (size_t j=1; j<size; ++j)
          if (j%Ny == Ny-1)
              Laplacian.coeffRef(j,j-1)           = 0.0;
}


void BaseLaplacian::SetRightSymmetry3()
{
  if (RightSymmetry == 1)
      for(size_t j=size-2*Ny; j<size-Ny; ++j)
          Laplacian.coeffRef(j+Ny,j)            += D1x;

  if (RightSymmetry == -1)
      for(size_t j=size-2*Ny; j<size-Ny; ++j)
          Laplacian.coeffRef(j+Ny,j)             = 0.0;

}

void BaseLaplacian::SetLeftSymmetry3()
{
  if (LeftSymmetry == 1)
      for(size_t j=0; j<Ny; ++j)
        Laplacian.coeffRef(j,j+Ny)              += D1x;

  if (LeftSymmetry == -1)
      for(size_t j=0; j<Ny; ++j)
          Laplacian.coeffRef(j,j+Ny)             = 0.0;
}

























void BaseLaplacian::SetRightSymmetry5()
{
  if (RightSymmetry == 1){
      for(size_t j=size-Ny; j<size; ++j){
          Laplacian.coeffRef(j,j-Ny)         += D1x ;
          Laplacian.coeffRef(j,j-2*Ny)       += D2x ;
          }
      for(size_t j=size-2*Ny; j<size-Ny; ++j)
          Laplacian.coeffRef(j,j+Ny)         += D2x ;
      }

  if (RightSymmetry == -1){
      for(size_t j=size-2*Ny; j<size-Ny; ++j){
          Laplacian.coeffRef(j,j-Ny)         = 0.0 ;
          Laplacian.coeffRef(j,j-2*Ny)       = 0.0 ;
          }
      for(size_t j=size-2*Ny; j<size-Ny; ++j)
          Laplacian.coeffRef(j,j+Ny)         = 0.0 ;
      }

}

void BaseLaplacian::SetLeftSymmetry5()
{
  std::cout<<"YOYOYOY\n";
  if (LeftSymmetry == 1){
      for(size_t j=0; j<Ny; ++j){
        Laplacian.coeffRef(j,j+Ny)           *= 2.0 ;
        Laplacian.coeffRef(j,j+2*Ny)         *= 2.0 ;
      }
      for(size_t j=Ny; j<2*Ny; ++j)
          Laplacian.coeffRef(j,j-Ny)         += D2x ;
  }


  if (LeftSymmetry == -1){
      for(size_t j=0; j<Ny; ++j){
        Laplacian.coeffRef(j,j+Ny)           = 0.0 ;
        Laplacian.coeffRef(j,j+2*Ny)         = 0.0 ;
      }
      for(size_t j=Ny; j<2*Ny; ++j)
          Laplacian.coeffRef(j,j-Ny)         = 0.0 ;
  }


}


void BaseLaplacian::SetTopSymmetry5()
{
  if (TopSymmetry == 1)
      for (size_t j=0; j<size; ++j)
          if (j%Ny == Ny-1){
              Laplacian.coeffRef(j,j-1)     += D1y ;
              Laplacian.coeffRef(j,j-2)     += D2y ;
              Laplacian.coeffRef(j-1,j)     += D2y ;
              }

  if (TopSymmetry == -1)
      for (size_t j=0; j<size; ++j)
          if (j%Ny == Ny-1){
              Laplacian.coeffRef(j,j-1)      = 0.0 ;
              Laplacian.coeffRef(j,j-2)      = 0.0 ;
              Laplacian.coeffRef(j-1,j)      = 0.0 ;
              }
}





void BaseLaplacian::SetBottomSymmetry5()
{
  if (BottomSymmetry == 1)
      for (size_t j=0; j<size; ++j)
          if (j%Ny == 0){
              Laplacian.coeffRef(j,j+1)     += D1y ;
              Laplacian.coeffRef(j,j+2)     += D2y ;
              Laplacian.coeffRef(j+1,j)     += D2y ;
              }

  if (BottomSymmetry == -1)
      for (size_t j=0; j<size; ++j)
          if (j%Ny == 0){
              Laplacian.coeffRef(j,j+1)      = 0.0 ;
              Laplacian.coeffRef(j,j+2)      = 0.0 ;
              Laplacian.coeffRef(j+1,j)      = 0.0 ;
              }
}
