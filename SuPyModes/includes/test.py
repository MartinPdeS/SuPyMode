from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
from SuPyModes.includes.EigenSolver  import EigenSolving
from SuPyModes.Geometry          import Geometry, Fused2, Circle
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

import time
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations


def PlotModes():
    for i in range(len(ITRList)):
        EigenVectors, EigenValues = A.GetSlice(i)

        fig, ax = plt.subplots(1,nMode)
        for j in range(nMode):
            array = EigenVectors[j,...]
            print(array.shape)
            temp = np.concatenate((array[:,::-1], array), axis=1)
            array = array.reshape([Nx,Ny])
            ax[j].pcolormesh(array)
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


nMode = 4

Nx = 50
Ny = 50

Clad = Fused2(Radius =  620.5, Fusion  = 1, Index   = 1)

Core0 = Circle( Position=(0,0), Radi = 15.2, Index = Fused_silica(1.55)+0.005 )
Xbound  = [-100, 100]
Ybound  = [-100, 100]
DeltaX = abs(Xbound[0]-Xbound[1])
DeltaY = abs(Ybound[0]-Ybound[1])

Geo = Geometry(Objects = [Clad],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)


ITRList = [1.0]

A = EigenSolving(Geo.mesh, Geo.mesh, nMode, 1550000, 1e-18, True)
A.dx = DeltaX/(Nx-1)

A.dy =  DeltaY/(Ny-1)
print(A.dx, A.dy)

A.Lambda          = 1.0

A.ComputeLaplacian()
A.TopSymmetry    = 0
A.BottomSymmetry = 0
A.LeftSymmetry   = 0
A.RightSymmetry  = 0

#A.PringLaplacian()

A.LoopOverITR(ITR = ITRList, alpha = -39.478, ExtrapolationOrder = 2)#-83.44

#A.SortModesFields()

def PlotIndex():
    index = A.ComputingIndices()

    for i in range(nMode):
        plt.plot(ITRList, index.reshape(len(ITRList),nMode)[:,i],'-', label=f'Mode: {i}')
    plt.grid()
    plt.legend()
    plt.show()





#PlotAdiabatic()
#PlotCoupling()
#PlotIndex()
PlotModes()




# -
