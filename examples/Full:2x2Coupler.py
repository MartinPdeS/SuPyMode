from SuPyModes.Geometry          import Geometry, Fused2, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(precision=3, linewidth=500)


nMode   = 6
Xbound  = [-100, 100]
Ybound  = [-100,100]
Nx      = 50
Ny      = 80

Clad = Fused2(Radius =  620.5, Fusion  = 0.95, Index   = Fused_silica(1.55))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )


Geo = Geometry(Objects = [Clad],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)

Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000, nMode=8, sMode=6)

#Sol.LeftSymmetry = 1
Sol.BottomSymmetry = 1

"""
Bottom -> Left
Left -> Bottom
Top -> Right

"""

SuperModes = Sol.GetModes(wavelength = 1.55*1.5,
                          Nstep      = 10,
                          ITRi       = 1,
                          ITRf       = 0.1)

SuperModes.Plot(Input=['Fields'], iter=[0])
#SuperModes.Plot(Input=['Fields'], iter=[0])
# SuperModes.Plot(Input=['Fields'], iter=[100])
# SuperModes.Plot(Input=['Fields'], iter=[150])
# SuperModes.Plot(Input=['Fields'], iter=[200])
# SuperModes.Plot(Input=['Fields'], iter=[250])
