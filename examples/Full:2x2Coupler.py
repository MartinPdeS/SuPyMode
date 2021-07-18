from SuPyModes.Geometry          import Geometry, Fused2, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(precision=3, linewidth=500)
Nx = 35
Ny = 30

Clad = Fused2(Radius =  620.5, Fusion  = 1, Index   = 1.0)

Core0 = Circle( Position=(0,0), Radi = 4.2, Index = Fused_silica(1.55)+0.005 )


Geo = Geometry(Objects = [Clad],
               Xbound  = [-100, 100],
               Ybound  = [-100, 100],
               Nx      = Nx,
               Ny      = Ny)

Sol = SuPySolver(Coupler=Geo)




SuperModes = Sol.GetModes(wavelength = 1.0,
                          Nstep      = 1,
                          Nsol       = 7,
                          debug      = False,
                          ITRi       = 1,
                          ITRf       = 0.99,
                          tolerance  = 1e-20,
                          error      = 2,
                          Xsym       = 1,
                          Ysym       = 0 )


SuperModes.Plot(Input=['Fields'])
