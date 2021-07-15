from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
from SuPyModes.utils             import *
from SuPyModes.perso             import FiberA, FiberB

A, B = FiberA(wavelength=1.3), FiberB(wavelength=1.3)

Capillary = Circle( Position = [0,0], Radi = 140, Index = 1.433,  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55), debug='WARNING')


Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )


Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


Geo = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-150, 150],
               Ybound  = [-150, 150],
               Nx      = 80,
               Ny      = 80,
               debug   = 'INFO',
               Length  = None)


Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength = 1.3,
                          Nstep      = 100,
                          Nsol       = 5,
                          ITRi       = 1,
                          ITRf       = 0.05,
                          tolerance  = 1e-20,
                          error      = 3,
                          Xsym       = 0,
                          Ysym       = 0 )


SuperModes.Plot(Input=['Index','Fields'], iter=[0], nMax=4)
