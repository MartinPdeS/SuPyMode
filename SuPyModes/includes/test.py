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

def PlotModes(iter):
    for i in iter:
        EigenVectors, EigenValues = A.GetSlice(i)

        fig, ax = plt.subplots(1,nMode)
        for j in range(nMode):
            array = EigenVectors[j,...]

            ax[j].pcolormesh(array.T, cmap='seismic')
            ax[j].set_aspect('equal')


            ax[j].set_title(f"beta: { EigenValues[j]:.5f} \n norm: {np.sum(np.abs(array)**2)}", fontsize=8)

        plt.tight_layout()
        plt.draw()

def PlotCoupling():

    Coupling = A.ComputingCoupling()

    a = list(combinations(range(nMode), 2))

    fig, ax = plt.subplots(1,1)
    for i,j in a:
        ax.plot(ITRList[1:], np.abs(Coupling[:,j,i]), label=f'{i}/{j}')

    plt.legend()
    plt.grid()
    plt.draw()


def PlotAdiabatic():

    Adiabatic = A.ComputingAdiabatic()

    print(Adiabatic)

    a = list(combinations(range(nMode), 2))

    fig, ax = plt.subplots(1,1)
    for i,j in a:
        ax.semilogy(ITRList[1:], Adiabatic[:,j,i], label=f'{i}/{j}')

    plt.legend()
    plt.grid()
    plt.draw()


def PlotIndex():
    index = A.ComputingIndices()

    for i in range(nMode):
        plt.plot(ITRList, index.reshape(len(ITRList),nMode)[:,i],'-', label=f'Mode: {i}')
    plt.grid()
    plt.legend()
    plt.draw()



nMode   = 5
Xbound  = [-100, 100]
Ybound  = [-100, 100]
Nx      = 90
Ny      = 90


Clad = Fused2(Radius =  62.5, Fusion  = 0.95, Index   = Fused_silica(1.55))

Core0 = Circle( Position=Clad.C[0], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position=Clad.C[1], Radi = 4.2, Index = Fused_silica(1.55)+0.005 )


Geo = Geometry(Objects = [Clad, Core0, Core1],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)

ITRList = np.linspace(1,0.9,10)

A = EigenSolving(Geo.mesh, Geo.Gradient().T.ravel(), nMode+3, nMode, 10000, 1e-18)

A.dx = Geo.Axes.dx

A.dy =  Geo.Axes.dy

A.Lambda          = 1.55*1.5

A.ComputeLaplacian()
A.TopSymmetry    = 0
A.BottomSymmetry = 0
A.LeftSymmetry   = 1
A.RightSymmetry  = 0

A.LoopOverITR(ITR = ITRList, ExtrapolationOrder = 3)

A.SortModesFields()


#PlotIndex()
PlotModes(iter=[0])
PlotCoupling()
PlotAdiabatic()




plt.show()


# -
