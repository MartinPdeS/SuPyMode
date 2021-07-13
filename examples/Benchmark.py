from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
import time

unit = 1e-0

Clad = Circle(Radi =  62.5*unit, Position = (0,0), Index = Fused_silica(1.0))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2*unit, Index = Fused_silica(1.0)+0.005 )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [-70*unit, 70*unit],
               Ybound  = [-70*unit, 70*unit],
               Nx      = 50,
               Ny      = 50)


Sol = SuPySolver(Coupler=Geo)


start = time.time()


SuperModes = Sol.GetModes(wavelength = 1*unit,
                          Nstep      = 1,
                          Nsol       = 5,
                          debug      = False,
                          ITRf       = 0.9,
                          Xsym       = 0,
                          Ysym       = 0 )

end = time.time()
print('time: ', end - start)

SuperModes.Plot(Input = ['Fields'], iter=[0])
