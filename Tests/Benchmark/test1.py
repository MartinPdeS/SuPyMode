from SuPyModes.Geometry          import Geometry, Circle, Fused2
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

import numpy as np
import matplotlib.pyplot as plt
"""
FIGURE 2.5 SBB_____________________________________________________

"""

Clad = Circle( Position = (0,0),
               Radi     = 62.5,
               Index    = Fused_silica(1.31)
               )

Core0 = Circle( Position = (0,0),
                Radi     = 4.1,
                Index    = Fused_silica(1.31)+0.005
                )

SMF28 = Geometry(Objects = [Clad, Core0],
                 Xbound  = [-90, 90],
                 Ybound  = [-90, 90],
                 Nx      = 50,
                 Ny      = 50,
                 GConv   = 0)

#SMF28.CreateMesh()
SMF28.Plot()

sMode=5

Sol = SuPySolver(Coupler    = SMF28,
                 Tolerance  = 1e-10,
                 MaxIter    = 10000,
                 nMode      = 2*sMode+1,
                 sMode      = sMode,
                 Error      = 2)

ITRList = np.linspace(1,0.1,300)
#ITRList = np.ones(30)
Sol.CppSolver.Lambda = 1.55
Sol.CppSolver.dx = SMF28.Axes.dx
Sol.CppSolver.dy = SMF28.Axes.dy
Sol.CppSolver.ComputeLaplacian(2)


Sol.CppSolver.LoopOverITR(ITRList, 3)

Field0, beta0 = Sol.CppSolver.GetSlice(199)

print(beta0)


fig, axes = plt.subplots(1,sMode)
# if not isinstance(axes, list):
#     axes = [axes]

for i in range(sMode):
    axes[i].pcolormesh(Field0[i])

plt.show()


"""
Note that the bad correspondence are due to the fact that we do not break the symmetries and
as such the modes are rotating from one slice to the other!

"""
