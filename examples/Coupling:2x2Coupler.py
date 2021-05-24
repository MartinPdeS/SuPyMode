from SuPyModes.Geometry          import Geometry, Fused2, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

Clad = Fused2(Radius =  62.5,
              Fusion  = 1,
              Index   = Fused_silica(1.55))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Geo = Geometry(Objects = [Clad, Core0, Core1],
               Xbound  = [0, 120],
               Ybound  = [0, 120],
               Nx      = 50,
               Ny      = 50,
               Xsym    = 1,
               Ysym    = 1)

Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength=1.55, Nstep=100, Nsol=6, debug=False, ITRf=0.02, tolerance=1e-20, error=3 )

SuperModes.Plot(['Adiabatic', 'Fields'], nMax = 4, iter=-1)
