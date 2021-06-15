from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica


Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))


Clad0 = Circle( Position=Clad.C[0], Radi = 12, Index = Fused_silica(1.55)+0.003 )

Clad1 = Circle( Position=Clad.C[1], Radi = 12, Index = Fused_silica(1.55)+0.003 )

Clad2 = Circle( Position=Clad.C[2], Radi = 10, Index = Fused_silica(1.55)+0.003 )


Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core2 = Circle( Position=Clad.C[2], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Geo = Geometry(Objects = [Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-120, 120],
               Ybound  = [-110, 130],
               Nx      = 10,
               Ny      = 10)

#Geo.Plot()

Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 100,
                          Nsol       = 7,
                          debug      = False,
                          ITRi       = 1,
                          ITRf       = 0.05,
                          tolerance  = 1e-20,
                          error      = 3,
                          Xsym       = 0,
                          Ysym       = 0 )

SuperModes.Plot(Input=['All'], nMax=3)
"""
SuperModes.Save(Directory  = 'lol.pdf',
                Input      = ['Index', 'Coupling', 'Adiabatic', 'Fields'],
                nMax       = 4)
"""
