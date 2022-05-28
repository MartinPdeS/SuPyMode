from SuPyMode.Geometry          import Geometry, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from SuPyMode.fibers            import *
from PyOptik                    import ExpData
import numpy as np


FiberA = Fiber_DCF1300S_20(Wavelength=1.55)
FiberB = Fiber_DCF1300S_33(Wavelength=1.55)

nMode   = 6
Xbound  = [-150, 150]
Ybound  = [-150, 150]
Nx      = 120
Ny      = 120

Index = ExpData('FusedSilica').GetRI(1.55e-6)

Clad = Fused4(Radius =  62.5, Fusion  = 0.95, Index   = Index)


Core0 = FiberA.Get(Clad.C[0])
Core1 = FiberB.Get(Clad.C[1])
Core2 = FiberA.Get(Clad.C[2])
Core3 = FiberB.Get(Clad.C[3])

Geo = Geometry(Clad    = Clad,
               Objects = [*Core0, *Core1, *Core2, *Core3],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)

#Geo.Rotate(45)

Geo.Plot()
dsa

Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000, nMode=8, sMode=5)

SuperSet = Sol.GetModes(wavelength        = 1.55,
                          Nstep           = 300,
                          ITRi            = 1,
                          ITRf            = 0.05,
                          Sorting         = 'Index',
                          RightSymmetry   = 0,
                          LeftSymmetry    = 0,
                          TopSymmetry     = 0,
                          BottomSymmetry  = 0
                          )

#SuperSet.PlotFields(iter=-1)
#SuperSet.Plot(Input=['Index', 'Adiabatic'], iter=[-1])
SuperSet.ExportPDF(Directory='ABAB', iter=[0, 200, 290])
