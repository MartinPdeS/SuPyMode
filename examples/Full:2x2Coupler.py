from SuPyModes.Geometry          import Geometry, Fused2, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(precision=3, linewidth=500)
Nx = 80
Ny = 80

Clad = Fused2(Radius =  62.5, Fusion  = 1, Index   = Fused_silica(1.55))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )


Geo = Geometry(Objects = [Clad, Core0, Core1],
               Xbound  = [-100, 100],
               Ybound  = [-100, 100],
               Nx      = Nx,
               Ny      = Ny)

Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 300,
                          Nsol       = 5,
                          debug      = False,
                          ITRi       = 1,
                          ITRf       = 0.1,
                          tolerance  = 1e-20,
                          error      = 2,
                          Xsym       = 0,
                          Ysym       = 0 )

SuperModes.Plot(Input=['Adiabatic', 'Index'])
# SuperModes.Plot(Input=['Fields'], iter=[50])
# SuperModes.Plot(Input=['Fields'], iter=[100])
# SuperModes.Plot(Input=['Fields'], iter=[150])
# SuperModes.Plot(Input=['Fields'], iter=[200])
# SuperModes.Plot(Input=['Fields'], iter=[250])
