from SuPyMode.Geometry          import Geometry, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica


Wavelength = 1.55e-6
Index = FusedSilica.GetRI(Wavelength)

Clad = Circle( Position=(0,0), Radius = 62.5, Index = Index)

Core = Circle( Position=Clad.C[0], Radius = 4.1, Index = Index+0.005 )

Geo = Geometry(Clad    = Clad,
               Objects = [Core],
               Xbound  = [-70, 70],
               Ybound  = [-70, 70],
               Nx      = 50,
               Ny      = 50)

#Geo.Plot()

Sol = SuPySolver(Coupler=Geo, Tolerance=1e-8, MaxIter = 10000)

Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)


Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': 0, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nMode           = 6,
             sMode           = 4 )

Set = Sol.GetSet()

Set.PlotFields([0, -1])

Set.PlotAdiabatic()
