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
                 Nx      = 80,
                 Ny      = 80,
                 GConv   = 0)


sMode=8

Sol = SuPySolver(Coupler    = SMF28,
                 Tolerance  = 1e-8,
                 MaxIter    = 1000,
                 nMode      = sMode,
                 sMode      = sMode,
                 Error      = 2)

ITRList = np.linspace(1,0.99,1)

Sol.CppSolver.Lambda = 1.55
Sol.CppSolver.dx = SMF28.Axes.dx
Sol.CppSolver.dy = SMF28.Axes.dy
Sol.CppSolver.ComputeLaplacian(2)


Sol.CppSolver.LoopOverITR_(ITRList, 1)

Field0, beta0 = Sol.CppSolver.GetSlice(0)

print(Field0.shape)


fig, axes = plt.subplots(1,sMode)
for i in range(sMode):
    axes[i].pcolormesh(Field0[i])

plt.show()


"""
Note that the bad correspondence are due to the fact that we do not break the symmetries and
as such the modes are rotating from one slice to the other!

"""
