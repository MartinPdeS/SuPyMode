vector<float> SD1A2 = {                           -0.5, +0.0, +0.5                                   };
vector<float> SD1A4 = {            +1.0/12.0, -2.0/3.0, +0.0, +2.0/3.0, -1.0/12.0                    };
vector<float> SD1A6 = { -1.0/60.0, +3.0/20.0, -3.0/4.0, +0.0, +3.0/4.0, -3.0/20.0 +1.0/60.0          };


vector<float> SD2A2 = {                           +1.0, -2.0,       +1.0                             };
vector<float> SD2A4 = {            -1.0/12.0, +4.0/3.0, -2.5,       +4.0/3.0, -1.0/12.0              };
vector<float> SD2A6 = { +1.0/90.0, -3.0/20.0, +3.0/2.0, -49.0/18.0, +3.0/2.0, -3.0/20.0 +1.0/60.0    };


vector<float> FD2A1 = { +1.0,      -2.0,       +1.0                                                  };
vector<float> FD2A2 = { +2.0,      -5.0,       +4.0,     -1.0                                        };
vector<float> FD2A3 = { 35.0/12.0, -26.0/3.0, +19.0/2.0, -14.0/3.0, 11.0/12.0                        };


//
// MSparse
// EigenSolving::Points3Laplacian_(){
//
//     size_t size = Nx*Ny;
//
//     MSparse Laplacian = MSparse(size,size);
//
//     for(uint j=0; j<size; ++j)
//         Laplacian.insert(j,j) = SD2A2[1]*( 1./pow(dx,2) + 1./pow(dx,2) );
//
//     for(uint j=0; j<size-1; ++j){
//       Laplacian.insert(j+1,j) = SD2A2[0]/pow(dx,2);
//       Laplacian.insert(j,j+1) = SD2A2[2]/pow(dx,2);
//       }
//
//     for(uint j=0; j<size-Nx; ++j){
//         Laplacian.insert(j+Nx,j) = SD2A2[0]/pow(dy,2);
//         Laplacian.insert(j,j+Nx) = SD2A2[2]/pow(dy,2);
//         }
//
//     for(uint j=1; j<size-1; ++j){
//         if (j%Nx == 0)   {Laplacian.coeffRef(j-1,j) = 0.0;}
//         if (j%Nx == Nx-1){Laplacian.coeffRef(j+1,j) = 0.0;}
//
//     }
//
//     return Laplacian;
// }
