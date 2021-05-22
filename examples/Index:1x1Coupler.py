from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

Clad = Circle(Radi =  62.5, Position = (0,0), Index = Fused_silica(1.0))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.0)+0.005 )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [0, 120],
               Ybound  = [0, 120],
               Nx      = 100,
               Ny      = 100,
               Xsym    = 1,
               Ysym    = 1)

Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength=1.0, Nstep=200, Nsol=8, debug=True )

SuperModes.Plot('Index')
