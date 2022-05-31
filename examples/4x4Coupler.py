from SuPyMode.Geometry          import Geometry, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from PyOptik                    import ExpData
import numpy as np


nMode   = 6
Xbound  = [-150, 0]
Ybound  = [-150, 0]
Nx      = 30
Ny      = 30

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


Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000)

SuperSet = Sol.GetSuperSet(Wavelength      = 1.55,
                          Nstep           = 300,
                          ITRi            = 1,
                          ITRf            = 0.05,
                          Sorting         = 'Index',
                          Symmetries      = {'Right': -1, 'Left': 0, 'Top': -1, 'Bottom': 0},
                          nMode=8,
                          sMode=5
                          )


#Sol.GetCoupling_()
Scene0 = SuperSet.PlotFields(iter=-1)

#SuperSet.PlotPropagation(Modes = [0,1])

#SuperSet.Plot(Input=['Coupling'])
#SuperSet.ExportPDF(Directory='4x4_SMF28_Hybrid_Ax_Ay', iter=[0, 100, 200, 290])
