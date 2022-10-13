"""
1x1 Coupler
===========

.. image:: ../../images/Example1/Geometry.png
   :width: 400
   :align: center
"""

# sphinx_gallery_thumbnail_path = '../images/Example1/Geometry.png'


from FiberFusing                import Geometry, Circle, BackGround
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica

Wavelength = 1.55e-6
Index = FusedSilica.GetRI(Wavelength)

Air = BackGround(Index=1) 

Clad = Circle( Position=(0,0), Radius = 62.5, Index = Index)

Core = Circle( Position=Clad.C[0], Radius = 4.1, Index = Index+0.005 )

Geo = Geometry(Objects = [Air, Clad, Core],
               Xbound  = [-70, 70],
               Ybound  = [-70, 70],
               Nx      = 180,
               Ny      = 180)

Geo.Plot().Show()

Sol = SuPySolver(Geometry=Geo, Tolerance=1e-8, MaxIter = 10000)

Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)


Sol.AddModes(Sorting         = 'Field',
             Symmetries      = {'Right': 0, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nComputedMode   = 8,
             nSortedMode     = 5 )


Set = Sol.GetSet()

Set.Plot(Type='field', ITR=[1, 0.3]).Show()

Set.Plot(Type='index').Show()

Set.Plot(Type='adiabatic').Show()