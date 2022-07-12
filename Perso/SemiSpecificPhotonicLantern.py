import numpy as np
from SuPyMode.Geometry          import Geometry, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica
from SuPyMode.fibers            import *

Wavelength = 1.55e-6
Index = FusedSilica.GetRI(Wavelength)


Angle = [0, 120, 240, 300, 360]


FiberA = Fiber_DCF1300S_20(Wavelength=1.55)
FiberB = Fiber_DCF1300S_33(Wavelength=1.55)
FiberC = Fiber_New(Wavelength=1.55)
#FiberA, FiberB = FiberB, FiberA


Clad = Circle(Position=(0,0), Radius = 100, Index = Index)

Cores = []

Cores += FiberA.Get(Position=(0,0))

for angle in Angle:
    P = (50*np.cos(angle*np.pi/180), 50*np.sin(angle*np.pi/180))
    Cores += FiberA.Get(Position=P)


for angle in [60]:
    P = (50*np.cos(angle*np.pi/180), 50*np.sin(angle*np.pi/180))
    Cores += FiberB.Get(Position=P)

for angle in [180]:
    P = (50*np.cos(angle*np.pi/180), 50*np.sin(angle*np.pi/180))
    Cores += FiberC.Get(Position=P)

Geo = Geometry(Clad    = Clad,
               Objects = Cores,
               Xbound  = [-120, 120],
               Ybound  = [-120, 120],
               Nx      = 120,
               Ny      = 120)

Geo.Rotate(90)

Geo.Plot()

Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000)

Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)


Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': 0, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nMode           = 10,
             sMode           = 8 )
#
# Sol.AddModes(Sorting         = 'Index',
#              Symmetries      = {'Right': -1, 'Left': 0, 'Top': 1, 'Bottom': 0},
#              nMode           = 7,
#              sMode           = 4 )
#
# Sol.AddModes(Sorting         = 'Index',
#              Symmetries      = {'Right': -1, 'Left': 0, 'Top': -1, 'Bottom': 0},
#              nMode           = 7,
#              sMode           = 4 )
#
#
# Sol.AddModes(Sorting         = 'Index',
#              Symmetries      = {'Right': 1, 'Left': 0, 'Top': -1, 'Bottom': 0},
#              nMode           = 7,
#              sMode           = 4 )


Set = Sol.GetSet()
#
Set.PlotFields(ITR=[0.9, 0.08])
#
# Set.PlotAdiabatic()



# Set.ExportPDF(Directory='TestNewNew', ITR=[0.9, 0.06])
