from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
import time

Clad = Circle(Radi =  62.5, Position = (0,0), Index = Fused_silica(1.0))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.0)+0.005 )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [-70, 70],
               Ybound  = [-70, 70],
               Nx      = 80,
               Ny      = 80)

#Geo.Plot()

Sol = SuPySolver(Coupler=Geo)


start = time.time()


SuperModes = Sol.GetModes(wavelength = 1.0,
                          Nstep      = 10,
                          Nsol       = 10,
                          debug      = False,
                          ITRf       = 0.9,
                          Xsym       = 0,
                          Ysym       = 0 )

end = time.time()
print('compute time:', end - start)

SuperModes.Plot(Input = ['All'], iter=[0])
