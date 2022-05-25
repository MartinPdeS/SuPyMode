from SuPyMode.Geometry          import Geometry, Fused3, Circle
from SuPyMode.Solver            import SuPySolver
from PyOptik                    import ExpData
import numpy as np


nMode   = 6
Xbound  = [-150, 150]
Ybound  = [-150, 150]
Nx      = 80
Ny      = 80

Wavelength = 1.55e-6
Index = ExpData('FusedSilica').GetRI(Wavelength)

Clad = Fused4(Radius =  62.5, Fusion  = 0.95, Index   = Index)

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Index+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Index+0.005 )

Core2 = Circle( Position=Clad.C[2], Radi = 4.2, Index = Index+0.005 )


Geo = Geometry(Objects = [Clad, Core0, Core1, Core2],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)

Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000, nMode=8, sMode=6)

SuperModes = Sol.GetModes(wavelength      = 1.55,
                          Nstep           = 10,
                          ITRi            = 1,
                          ITRf            = 0.1,
                          Sorting         = 'Field',
                          RightSymmetry   = 0,
                          LeftSymmetry    = 0,
                          TopSymmetry     = 0,
                          BottomSymmetry  = 0
                          )


SuperModes.Plot(Input=['Fields'], iter=[0])
