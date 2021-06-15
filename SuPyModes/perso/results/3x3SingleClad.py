from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica


Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Core0 = Circle( Position=Clad.C[0], Radi = 4.1, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.1, Index = Fused_silica(1.55)+0.005 )

Core2 = Circle( Position=Clad.C[2], Radi = 4.1, Index = Fused_silica(1.55)+0.005 )

Geo = Geometry(Objects = [Clad, Core0, Core1, Core2],
               Xbound  = [-120, 120],
               Ybound  = [-110, 130],
               Nx      = 80,
               Ny      = 80)


Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 200,
                          Nsol       = 7,
                          debug      = False,
                          ITRi       = 1,
                          ITRf       = 0.05,
                          tolerance  = 1e-20,
                          error      = 3,
                          Xsym       = 0,
                          Ysym       = 0 )

SuperModes.SaveFig(Directory  = 'perso/SingleClad',
                   Input      = ['All'],
                   nMax       = 4,
                   iter       = -1,)
