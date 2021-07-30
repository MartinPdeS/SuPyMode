from SuPyModes.Geometry          import Geometry, Circle, Fused2
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

from SuPyModes.perso.fibers import Fiber_SMF28

A = Fiber_SMF28(wavelength=1.55)

Clad = Circle( Position = (0,0), Radi = 62.5, Index = Fused_silica(1.55) )

Core0 = Circle( Position = (0,0), Radi = 4.1, Index = Fused_silica(1.55)+0.005 )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [-70, 0],
               Ybound  = [-70, 0],
               Nx      = 100,
               Ny      = 100)


#Geo.Plot()

Sol = SuPySolver(Coupler    = Geo,
                 Tolerance  = 1e-40,
                 MaxIter    = 1000,
                 nMode      = 9,
                 sMode      = 7,
                 )

SuperModesSym = Sol.GetModes(wavelength = 1.55,
                          Nstep         = 300,
                          ITRi          = 1,
                          ITRf          = 0.1,
                          RightSymmetry = 1,
                          TopSymmetry   = 1,
                          Sorting       = 'Fields')

SuperModesSym.Plot(Input   = ['Fields','Adiabatic'], iter=[0,-1], nMax=7)
