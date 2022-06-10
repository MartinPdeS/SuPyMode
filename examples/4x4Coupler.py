

from SuPyMode.Geometry          import Geometry, Fused2, Circle
from SuPyMode.Solver            import SuPySolver
from PyOptik                    import ExpData
import numpy as np
import matplotlib.pyplot as plt

#np.set_printoptions(edgeitems=30, linewidth=1000,
#    formatter=dict(float=lambda x: "%.3e" % x))



nMode   = 6

N = Nx = Ny  = 40


Index = ExpData('FusedSilica').GetRI(1.55e-6)

Clad = Fused2(Radius =  62.5, Fusion  = 0.95, Index   = Index)

Core0 = Circle( Position=Clad.C[0], Radius = 4.2, Index = Index+0.005 )
Core1 = Circle( Position=Clad.C[1], Radius = 8.2, Index = Index+0.005 )


Geo = Geometry(Clad    = Clad,
               Objects = [Core0, Core1],
               Xbound  = [-150, 0],
               Ybound  = [-150, 150],
               Nx      = Nx//2,
               Ny      = Ny,
               BackGroundIndex = 1.)

Geo.Rotate(0)

#Geo.Plot()


Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000, Debug=True)

Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)

Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': 1, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nMode=6,
             sMode=5 )


Set = Sol.GetSet()
Set.ExportPDF()
#Set.PlotFields([0, -1])
#Set[0].PlotFields([0, -1])





# -
