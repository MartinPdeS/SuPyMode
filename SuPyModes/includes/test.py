from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
from SuPyModes.includes.EigenSolver  import EigenSolving
import time
import numpy as np
import matplotlib.pyplot as plt

Nx = 80
Ny = 80

Clad = Circle(Radi =  62.5, Position = (0,0), Index = Fused_silica(1.0))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.0)+0.005 )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [-100, 100],
               Ybound  = [-100, 100],
               Nx      = Nx,
               Ny      = Ny)

ITRList = np.linspace(1.0,0.9, 10)

A = EigenSolving(Geo.mesh, 5, 500, 1e-10)
A.dx = A.dy = 140/(Nx-1)
A.Lambda = 1.0/ITRList[1]

A.LoopOverITR(ITRList, -83.44)
Coupling = A.ComputingOverlap()
print(Coupling)
# for i in range(len(ITRList)):
#     EigenVectors, EigenValues = A.GetSlice(i)
#
#     fig, ax = plt.subplots(1,5)
#     for j in range(5):
#         ax[j].imshow(EigenVectors[j,...])
#         print(EigenVectors)
#         ax[j].set_title(f"beta: { EigenValues[j]:.2f}", fontsize=8)
#
#     plt.tight_layout()
#     plt.show()





# print(ITRList[1])
#
# EigenVectors, EigenValues =  A.ComputeEigen(-71)
#
# print(EigenValues.shape)
# fig, ax = plt.subplots(1,5)
# for j in range(5):
#
#     ax[j].imshow(EigenVectors[j,...])
#     ax[j].set_title(f"beta: { EigenValues[j]:.2f}", fontsize=8)
#
# plt.tight_layout()
# plt.show()






        #
