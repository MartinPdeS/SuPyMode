from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
from SuPyModes.includes.EigenSolver  import EigenSolving
from SuPyModes.Geometry          import Geometry, Fused2, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
np.set_printoptions(threshold=sys.maxsize)

def PlotModes():
    for i in range(len(ITRList)):
        EigenVectors, EigenValues = A.GetSlice(i)

        fig, ax = plt.subplots(1,nMode)
        for j in range(nMode):
            array = EigenVectors[j,...]
            #print(array.shape)
            #temp = np.concatenate((array[:,::-1], array), axis=1)
            #array = array.reshape([Nx,Ny])
            ax[j].pcolormesh(array, cmap='seismic')
            ax[j].set_aspect('equal')


            ax[j].set_title(f"beta: { EigenValues[j]:.5f}", fontsize=8)

        plt.tight_layout()
        plt.show()

def PlotCoupling():

    Coupling = A.ComputingCoupling()

    print(Coupling.shape)

    a = list(combinations(range(nMode), 2))

    fig, ax = plt.subplots(1,1)
    for i,j in a:
        ax.plot(ITRList[1:], np.abs(Coupling[:,j,i]), label=f'{i}/{j}')

    plt.legend()
    plt.grid()
    plt.show()


def PlotAdiabatic():

    Adiabatic = A.ComputingAdiabatic()

    print(Adiabatic)

    a = list(combinations(range(nMode), 2))

    fig, ax = plt.subplots(1,1)
    for i,j in a:
        ax.semilogy(ITRList[1:], np.abs(Adiabatic[:,j,i]), label=f'{i}/{j}')

    plt.legend()
    plt.grid()
    plt.show()


def PlotIndex():
    index = A.ComputingIndices()

    for i in range(nMode):
        plt.plot(ITRList, index.reshape(len(ITRList),nMode)[:,i],'-', label=f'Mode: {i}')
    plt.grid()
    plt.legend()
    plt.show()



nMode = 6

Nx = 200
Ny = 200
Xbound  = [-100, 100]
Ybound  = [-100, 100]

class Gradient1:
    def __init__(self, center):
        self.center = center

    def evaluate(self, Xmesh, Ymesh):
        rho = np.sqrt( (Xmesh-self.center[0])**2 + (Ymesh-self.center[1])**2 )+10
        return 1/rho


Clad = Fused2(Radius =  62.5, Fusion  = 0.3, Index   = Fused_silica(1.55))

Clad.Gradient = Gradient1

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )


Geo = Geometry(Objects = [Clad],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)

Geo.Plot()


"""
def Gradient(Index, Xmesh, Ymesh, center, order=1)
    rho = sqrt( (Xmesh-center[0])**2 + (Ymesh-center[1])**2 )
    mesh = Index * rho




DeltaX = abs(Xbound[0]-Xbound[1])
DeltaY = abs(Ybound[0]-Ybound[1])

ITRList = np.linspace(1,0.3,300)

A = EigenSolving(Geo.mesh, Geo.Gradient().T, nMode, 10000, 1e-18)
A.dx = DeltaX/(Nx-1)

A.dy =  DeltaY/(Ny-1)
print(A.dx, A.dy)

A.Lambda          = 1.55

A.ComputeLaplacian()
A.TopSymmetry    = 0
A.BottomSymmetry = 0
A.LeftSymmetry   = 0
A.RightSymmetry  = 0

A.LoopOverITR(ITR = ITRList, ExtrapolationOrder = 2)

A.SortModesFields()


#PlotIndex()
#PlotAdiabatic()
#PlotCoupling()
#PlotModes()
"""



# -
