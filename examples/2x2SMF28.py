from SuPyMode.Geometry          import Geometry, Fused2, Circle
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica


Wavelength = 1.55e-6
Index = FusedSilica.GetRI(Wavelength)

Clad = Fused2(Radius = 62.5, Fusion = 1, Index = Index)

Core0 = Circle( Position=Clad.C[0], Radius = 4.1, Index = Index+0.005 )

Core1 = Circle( Position=Clad.C[1], Radius = 4.1, Index = Index+0.005 )

Geo = Geometry(Clad    = Clad,
               Objects = [Core0, Core1],
               Xbound  = [-100, 0],
               Ybound  = [-100, 100],
               Nx      = 40,
               Ny      = 80)
Geo.Rotate(90)

#Geo.Plot()

Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000)

Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)


Sol.AddModes(Sorting         = 'Field',
             Symmetries      = {'Right': 1, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nMode           = 8,
             sMode           = 5 )

# Sol.AddModes(Sorting         = 'Index',
#              Symmetries      = {'Right': -1, 'Left': 0, 'Top': 0, 'Bottom': 0},
#              nMode           = 4,
#              sMode           = 2 )

Set = Sol.GetSet()

# Set.PlotFields([0, -1])

Set.PlotCoupling()
