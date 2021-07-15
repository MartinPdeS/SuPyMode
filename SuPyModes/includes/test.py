from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
from SuPyModes.includes.EigenSolver  import EigenSolving
import time
import numpy as np
import matplotlib.pyplot as plt

Nx = 80
Ny = 80
nMode = 2

Clad = Circle(Radi =  62.5, Position = (0,0), Index = Fused_silica(1.0))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.0)+0.005 )

Geo = Geometry(Objects = [Clad, Core0],
               Xbound  = [-80, 80],
               Ybound  = [-80, 80],
               Nx      = Nx,
               Ny      = Ny)

ITRList = np.linspace(1.0,0.3, 200)

A = EigenSolving(Geo.mesh, nMode+3, 1500, 1e-8)
A.dx = A.dy = 140/(Nx-1)
A.Lambda = 1.0

A.LoopOverITR(ITRList, -83.44, 1)

A.SortModesFields(nMode)
#A.SortModesIndex()

index = A.ComputingIndices()
print(index.shape)

for i in range(nMode):
    plt.plot(ITRList, index.reshape(len(ITRList),nMode)[:,i],'*-', label=f'Mode: {i}')
plt.grid()
plt.legend()
plt.show()




# for i in range(len(ITRList)):
#     EigenVectors, EigenValues = A.GetSlice(i)
#     print(EigenValues)
#
#     fig, ax = plt.subplots(1,nMode)
#     for j in range(nMode):
#         ax[j].imshow(EigenVectors[j,...])
#
#         print(EigenValues)
#         ax[j].set_title(f"beta: { EigenValues[j]:.5f}", fontsize=8)
#
#     plt.tight_layout()
#     plt.show()
#




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
