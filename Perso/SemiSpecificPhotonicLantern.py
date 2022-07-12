import numpy as np
from SuPyMode.Geometry          import Geometry, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica
from SuPyMode.fibers            import *

Wavelength = 1.55e-6
Index = FusedSilica.GetRI(Wavelength)


AngleFiberA = [60, 120, 240, 300, 360]
AngleFiberB = [0, 180]
AngleFiberC = [None]
A, B, C = Fiber_DCF1300S_20, Fiber_DCF1300S_33, Fiber_New


Config = [(0, A, B, C),
          (1, A, C, B),
          (2, C, B, A),
          (3, C, A, B),
          (4, B, A, C),
          (5, B, C, A)]


confNumber, FiberA, FiberB, FiberC = Config[5]




Clad = Circle(Position=(0,0), Radius = 100, Index = Index)

Cores = []


for angle in AngleFiberA:
    if angle is None:
        P = (0,0)
    else:
        P = (50*np.cos(angle*np.pi/180), 50*np.sin(angle*np.pi/180))
    Cores += FiberA(Wavelength=1.55).Get(Position=P)

for angle in AngleFiberB:
    if angle is None:
        P = (0,0)
    else:
        P = (50*np.cos(angle*np.pi/180), 50*np.sin(angle*np.pi/180))
    Cores += FiberB(Wavelength=1.55).Get(Position=P)

for angle in AngleFiberC:
    if angle is None:
        P = (0,0)
    else:
        P = (50*np.cos(angle*np.pi/180), 50*np.sin(angle*np.pi/180))
    Cores += FiberC(Wavelength=1.55).Get(Position=P)

Geo = Geometry(Clad    = Clad,
               Objects = Cores,
               Xbound  = [-120, 120],
               Ybound  = [-120, 120],
               Nx      = 120,
               Ny      = 120)

# Geo.Rotate(0)

# Geo.Plot()

Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000)

Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)


Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': 0, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nMode           = 20,
             sMode           = 17 )

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
# Set.PlotFields(ITR=[0.9, 0.08])
#
# Set.PlotAdiabatic()



Set.ExportPDF(Directory=f'Perso/17ModesHybridMSPL_config_{confNumber}', ITR=[0.9, 0.06])
