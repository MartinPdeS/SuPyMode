from SuPyMode.Geometry          import Geometry, Fused2, Circle
from SuPyMode.Solver            import SuPySolver
from SuPyMode.sellmeier         import Fused_silica
import numpy as np


nMode   = 6
Xbound  = [-50, 50]
Ybound  = [0,100]
Nx      = 80
Ny      = 80

Clad = Fused2(Radius =  62.5, Fusion  = 0.95, Index   = Fused_silica(1.55))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )


Geo = Geometry(Objects = [Clad, Core0, Core1],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)

Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000, nMode=8, sMode=6)

SuperModes = Sol.GetModes(wavelength      = 1.55*1.5,
                          Nstep           = 30,
                          ITRi            = 1,
                          ITRf            = 0.8,
                          Sorting         = 'Fields',
                          RightSymmetry   = 1,
                          LeftSymmetry    = 1,
                          TopSymmetry     = 0,
                          BottomSymmetry  = 1
                          )


SuperModes.Plot(Input=['Fields'], iter=[0])
