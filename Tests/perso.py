from SuPyMode.Geometry          import Geometry, Fused2, Fused3, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from PyOptik                    import ExpData
import numpy as np
import matplotlib.pyplot as plt
from matplotlib          import colors
import matplotlib

from SuPyMode.fibers import *

Index = ExpData('FusedSilica').GetRI(1.55e-6)



Index = ExpData('FusedSilica').GetRI(1.55e-6)


Clad = Fused4(Radius =  62.5, Fusion  = 0.99, Index   = Index)

Core0 = Circle( Position=Clad.C[0], Radius = 4.1, Index = Index+0.05 )
Core1 = Circle( Position=Clad.C[1], Radius = 4.1, Index = Index+0.05 )
Core2 = Circle( Position=Clad.C[2], Radius = 4.1, Index = Index+0.05 )
Core3 = Circle( Position=Clad.C[3], Radius = 4.1, Index = Index+0.05 )

Geo0 = Geometry(Clad    = Clad,
               Objects = [Core0, Core1, Core2, Core3],
               Xbound  = [-150, 0],
               Ybound  = [-150, 0],
               Nx      = 10,
               Ny      = 10,
               BackGroundIndex = 1.0)

Geo0.Rotate(90)


Geo0.Plot()


Sol = SuPySolver(Coupler=Geo0, Tolerance=1e-8, MaxIter = 10000)


Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)

Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': -1, 'Left': 0, 'Top': -1, 'Bottom': 0},
             nMode           = 4,
             sMode           = 2 )

print("1 Finished")

Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': 1, 'Left': 0, 'Top': -1, 'Bottom': 0},
             nMode           = 4,
             sMode           = 2 )
print("2 Finished")

Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': -1, 'Left': 0, 'Top': 1, 'Bottom': 0},
             nMode           = 4,
             sMode           = 2 )
print("3 Finished")

Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': 1, 'Left': 0, 'Top': 1, 'Bottom': 0},
             nMode           = 4,
             sMode           = 2 )
print("4 Finished")


Set = Sol.GetSet()
Set.PlotFields([0, -1])

Set.Plot("Adiabatic")
