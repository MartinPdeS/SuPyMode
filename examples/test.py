from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

Clad = Circle(Radi =  62.5, Position = (1000,0), Index = 1.01)


Geo = Geometry(Objects = [Clad],
               Xbound  = [-50, 50],
               Ybound  = [-50, 50],
               Nx      = 10,
               Ny      = 10)


Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength = 1,
                          Nstep      = 1,
                          Nsol       = 5,
                          debug      = False,
                          Xsym       = 0,
                          Ysym       = 0 )

SuperModes.Plot(Input=['Fields'])
