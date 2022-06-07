

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
Core1 = Circle( Position=Clad.C[1], Radius = 4.2, Index = Index+0.005 )


Geo = Geometry(Clad    = Clad,
               Objects = [Core0, Core1],
               Xbound  = [-150, 0],
               Ybound  = [-150, 150],
               Nx      = Nx//2,
               Ny      = Ny,
               BackGroundIndex = 1.)

Geo.Rotate(90)

#Geo.Plot()


Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000)

SuperSet = Sol.GetSuperSet(Wavelength      = 1.55,
                           Nstep           = 300,
                           ITRi            = 1,
                           ITRf            = 0.05,
                           Sorting         = 'Index',
                           Symmetries      = {'Right': 1, 'Left': 0, 'Top': 0, 'Bottom': 0},
                           nMode=5,
                           sMode=3
                           )

#SuperSet.PlotFields([-1])

#Mode = SuperSet.GetSuperposition(Amplitudes=[1,0,0,1,0,0])

#Mode.CouplerLength = 500
#Mode.ComputeAmpltiudes(MaxStep=0.1)


#Mode.PlotField(Slice=0)


"""


Mode.PlotFields(Slice=[0])

print(dir(SuperSet[0].CppSolver.GetMode(0)))
dsa
Fields = SuperSet[0].CppSolver.GetFields()
print(Fields.shape)


import matplotlib.pyplot as plt
Fig, ax = plt.subplots(1, 5)
for i in range(5):
    ax[i].pcolormesh(Fields[i,0,:,:])

plt.show()
#a = SuperSet[0].CppSolver.GetMode(2)

#SuperSet.PlotFields()
#Set.PlotPropagation(Modes = [0,1])

#SuperSet.Plot(Input=['Coupling'])
#SuperSet.ExportPDF(Directory='4x4_SMF28_Hybrid_Ax_Ay', iter=[0, 100, 200, 290])
"""
