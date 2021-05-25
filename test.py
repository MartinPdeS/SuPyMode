
from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

Clad = Circle(Radi =  62.5, Position = (0,0), Index = Fused_silica(1.55))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [0, 120],
               Ybound  = [0, 120],
               Nx      = 50,
               Ny      = 50)
Geo.Plot()
Sol = SuPySolver(Coupler=Geo)

SuperSet = Sol.GetModes(wavelength=1.55, Nstep=15, Nsol=7, debug=False, ITRi=1, ITRf=0.05, tolerance=1e-16, Xsym = 1, Ysym    = 1)
SuperSet.Plot('Adiabatic', iter=34, nMax=5)
