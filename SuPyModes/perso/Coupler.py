A, B = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = Fused_silica(1.55)  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


BBA_silica = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-150, 150],
               Ybound  = [-150, 150],
               Nx      = 150,
               Ny      = 150)



B, A = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = Fused_silica(1.55)  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


AAB_silica = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-150, 150],
               Ybound  = [-150, 150],
               Nx      = 150,
               Ny      = 150)




A, B = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = 1.433  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


BBA_fluoride = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-150, 150],
               Ybound  = [-150, 150],
               Nx      = 150,
               Ny      = 150)


B, A = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = 1.433  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


AAB_fluoride = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
                       Xbound  = [-150, 150],
                       Ybound  = [-150, 150],
                       Nx      = 150,
                       Ny      = 150)



B, A = Fiber_SMF28(wavelength=1.55), Fiber_SMF28(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = 1.433  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


SMF_fluoride = Geometry(Objects = [Capillary, Clad, Core0, Core1, Core2],
                       Xbound  = [-150, 150],
                       Ybound  = [-150, 150],
                       Nx      = 150,
                       Ny      = 150)



B, A = Fiber_SMF28(wavelength=1.55), Fiber_SMF28(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = Fused_silica(1.55)  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


SMF_silica = Geometry(Objects = [Capillary, Clad, Core0, Core1, Core2],
                       Xbound  = [-150, 150],
                       Ybound  = [-150, 150],
                       Nx      = 150,
                       Ny      = 150)
