from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

from SuPyModes.perso.fibers import Fiber_SMF28

A = Fiber_SMF28(wavelength=1.31)


Clad = Circle( Position = [0,0], Radi = 62.5, Index = A.nClad )

Core0 = Circle( Position = [0,0], Radi = A.rCore, Index = A.nCore )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [-80, 80],
               Ybound  = [-80, 80],
               Nx      = 10,
               Ny      = 10)


#Geo.Plot()

Sol = SuPySolver(Coupler    = Geo,
                 Tolerance  = 1e-20,
                 MaxIter    = 1000,
                 nMode      = 5,
                 sMode      = 3)

SuperModes = Sol.GetModes(wavelength = 1.31,
                          Nstep      = 1,
                          ITRi       = 1,
                          ITRf       = 0.99)

SuperModes.Plot(Input      = ['Fields'] )
