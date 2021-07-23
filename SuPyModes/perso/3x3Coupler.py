from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

from SuPyModes.perso.fibers import Fiber_DCF1300S_33, Fiber_DCF1300S_20, Fiber_SMF28

A, B = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = Fused_silica(1.55)  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )


Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


Geo = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-150, 150],
               Ybound  = [-150, 150],
               Nx      = 150,
               Ny      = 150)

Geo.Rotate(29)
Geo.Plot()

Sol = SuPySolver(Coupler    = Geo,
                 Tolerance  = 1e-20,
                 MaxIter    = 1000,
                 nMode      = 7,
                 sMode      = 5)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 300,
                          ITRi       = 1,
                          ITRf       = 0.05)


SuperModes.SaveFig(Directory  = 'test',
                   Input      = ['All'],
                   nMax       = 5)
