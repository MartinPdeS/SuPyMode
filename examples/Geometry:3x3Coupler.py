from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

Clad = Fused3(Radius =  62.5,
              Fusion  = 0.8,
              Index   = Fused_silica(1.55))

Clad.Plot()

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core2 = Circle( Position=Clad.C[2], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )


Geo = Geometry(Objects = [Clad, Core0, Core1, Core2],
               Xbound  = [0, 120],
               Ybound  = [-110, 130],
               Nx      = 100,
               Ny      = 100,
               Xsym    = 1,
               Ysym    = 0)

Geo.Plot()

Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength=1.55, Nstep=2, Nsol=3 )



#SuperModes.Plot('Adiabatic')

#SuperModes.Plot('Fields')

#SuperModes.Plot('Coupling',iter=0)
