#!/usr/bin/env python3

from SuPyModes.Geometry import Geometry, Circle
from SuPyModes.Solver import SuPySolver
from SuPyModes.sellmeier import Fused_silica
import time

Clad = Circle(Radi=62.5, Position=(0, 0), Index=Fused_silica(1.55))

Core0 = Circle(Position=Clad.C[0], Radi=4.2, Index=Fused_silica(1.55)+0.005)

Geo = Geometry(Objects=[Clad, Core0],
               Xbound=[-70, 70],
               Ybound=[-70, 70],
               Nx=150,
               Ny=150)

#Geo.Plot()

Sol = SuPySolver(Coupler=Geo,
                 Tolerance=1e-8,
                 MaxIter=10000,
                 nMode=8,
                 sMode=6)


start = time.time()


SuperModes = Sol.GetModes(wavelength=1.55,
                          Nstep=2,
                          ITRf=0.9)

end = time.time()
print('compute time:', end - start)

SuperModes.Plot(Input=['All'], iter=[0])
