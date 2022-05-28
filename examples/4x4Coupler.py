from SuPyMode.Geometry          import Geometry, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from PyOptik                    import ExpData
import numpy as np


nMode   = 6
Xbound  = [-150, 0]
Ybound  = [-150, 0]
Nx      = 80
Ny      = 80

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

Geo.Rotate(45)

#Geo.Plot()


Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000, nMode=8, sMode=5)

SuperSet = Sol.GetModes(wavelength      = 1.55,
                          Nstep           = 3,
                          ITRi            = 1,
                          ITRf            = 0.05,
                          Sorting         = 'Index',
                          RightSymmetry   = -1,
                          LeftSymmetry    = 0,
                          TopSymmetry     = -1,
                          BottomSymmetry  = 0
                          )

SuperSet.PlotFields(iter=-1)
#SuperSet.Plot(Input=['Index', 'Adiabatic'], iter=[-1])
#SuperSet.ExportPDF(Directory='4x4_SMF28_Hybrid_Ax_Ay', iter=[0, 100, 200, 290])
