from SuPyMode.Geometry          import Geometry, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from PyOptik                    import ExpData
import numpy as np


nMode   = 6
Xbound  = [-150, 150]
Ybound  = [-150, 150]
Nx      = 120
Ny      = 120

Index = ExpData('FusedSilica').GetRI(1.55e-6)

Clad = Fused4(Radius =  62.5, Fusion  = 0.95, Index   = Index)

Core0 = Circle( Position=Clad.C[0], Radius = 4.2, Index = Index+0.005 )

Core1 = Circle( Position=Clad.C[1], Radius = 4.2, Index = Index+0.005 )

Core2 = Circle( Position=Clad.C[2], Radius = 4.2, Index = Index+0.005 )

Core3 = Circle( Position=Clad.C[3], Radius = 4.2, Index = Index+0.005 )


Geo = Geometry(Clad    = Clad,
               Objects = [Core0, Core1, Core2, Core3],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)


#Geo.Plot()


Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000, nMode=8, sMode=5)

SuperModes = Sol.GetModes(wavelength      = 1.55,
                          Nstep           = 300,
                          ITRi            = 1,
                          ITRf            = 0.05,
                          Sorting         = 'Index',
                          RightSymmetry   = 0,
                          LeftSymmetry    = 0,
                          TopSymmetry     = 0,
                          BottomSymmetry  = 0
                          )

#SuperModes.PlotFields(iter=0)
#SuperModes.Plot(Input=['Index', 'Adiabatic'], iter=[-1])
SuperModes.ExportPDF(Directory='RodrigoMexicano', iter=0)
