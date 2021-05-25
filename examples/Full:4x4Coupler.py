from SuPyModes.Geometry          import Geometry, Circle, Fused4
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

Clad = Fused4(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core2 = Circle( Position=Clad.C[2], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core3 = Circle( Position=Clad.C[3], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Geo = Geometry(Objects = [Clad, Core0, Core1, Core2, Core3],
               Xbound  = [0, 160],
               Ybound  = [0, 160],
               Nx      = 70,
               Ny      = 70)

Geo.Plot()

Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 250,
                          Nsol       = 7,
                          debug      = False,
                          ITRi       = 1,
                          ITRf       = 0.05,
                          tolerance  = 1e-20,
                          error      = 3,
                          Xsym       = 1,
                          Ysym       = -1 )

SuperModes.Plot(Input = ['Index', 'Coupling', 'Adiabatic', 'Fields'], nMax=4)
