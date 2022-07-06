from SuPyMode.Geometry          import Geometry, Fused3, Circle
from SuPyMode.Solver            import SuPySolver
from SuPyMode.fibers            import *
from PyOptik                    import ExpData
import numpy as np


FiberA = Fiber_DCF1300S_20(Wavelength=1.55)
FiberB = Fiber_DCF1300S_33(Wavelength=1.55)
FiberC = Fiber_New(Wavelength=1.55)

Xbound  = [-130, 130]
Ybound  = [-130, 130]
Nx      = 130
Ny      = 130

Index = ExpData('FusedSilica').GetRI(1.55e-6)

Clad = Fused3(Radius = 62.5, Fusion = 0.95, Index   = Index)


Core0 = FiberA.Get(Clad.C[0])
Core1 = FiberB.Get(Clad.C[1])
Core2 = FiberC.Get(Clad.C[2])


Geo = Geometry(Clad    = Clad,
               Objects = [*Core0, *Core1, *Core2],
               Xbound  = Xbound,
               Ybound  = Ybound,
               Nx      = Nx,
               Ny      = Ny)

Geo.Rotate(45)

Geo.Plot()


Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000)
Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)

Sol.AddModes(Sorting    = 'Field',
             Symmetries = {'Right': 0, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nMode      = 8,
             sMode      = 4 )


Set = Sol.GetSet()

Set.PlotAdiabatic()

# SuperSet.ExportPDF(Directory='ABAB_xS_yA', iter=[0, 200, 290])
