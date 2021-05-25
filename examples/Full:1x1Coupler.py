from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

Clad = Circle(Radi =  62.5, Position = (0,0), Index = Fused_silica(1.0))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.0)+0.005 )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [0, 70],
               Ybound  = [0, 70],
               Nx      = 80,
               Ny      = 80)

Geo.Plot()

Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength = 1.0,
                          Nstep      = 250,
                          Nsol       = 10,
                          debug      = False,
                          Xsym       = 1,
                          Ysym       = -1 )

SuperModes.Plot(Input = ['Index', 'Coupling', 'Adiabatic', 'Fields'], nMax=4)
