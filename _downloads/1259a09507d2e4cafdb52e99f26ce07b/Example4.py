"""
4x4 Coupler
===========

.. image:: ../../images/Example4/Geometry.png
   :width: 400
   :align: center

.. image:: ../../images/Example4/Fields.png
   :width: 800
   :align: center


.. image:: ../../images/Example4/Index.png
   :width: 600
   :align: center


.. image:: ../../images/Example4/Adiabatic.png
   :width: 600
   :align: center
"""

# sphinx_gallery_thumbnail_path = '../images/Example4/Geometry.png'


from FiberFusing                import Geometry, Fused4, Circle, BackGround
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica

Air = BackGround(Index=1) 

Wavelength = 1.55e-6
Index = FusedSilica.GetRI(Wavelength)

Clad = Fused4(FiberRadius = 62.5, Fusion = 0.5, Index = Index)

Cores =  [ Circle(Position=Core, Radius=4.1, Index=Index+0.005) for Core in Clad.Cores]

Geo = Geometry(Objects = [Air, Clad] + Cores,
               Xbound  = [-150, 0],
               Ybound  = [-150, 0],
               Nx      = 80,
               Ny      = 80)

Geo.Plot().Show()

Sol = SuPySolver(Geometry=Geo, Tolerance=1e-8, MaxIter=10000)

Sol.CreateSuperSet(Wavelength=1.55, NStep=500, ITRi=1, ITRf=0.05)


Sol.AddModes(Sorting         = 'Field',
             Symmetries      = {'Right': 1, 'Left': 0, 'Top': 1, 'Bottom': 0},
             nComputedMode   = 3,
             nSortedMode     = 2 )

Sol.AddModes(Sorting         = 'Field',
             Symmetries      = {'Right': -1, 'Left': 0, 'Top': 1, 'Bottom': 0},
             nComputedMode   = 3,
             nSortedMode     = 2 )

Set = Sol.GetSet()

Set.Plot(Type='field', ITR=[1, 0.3]).Show()

Set.Plot(Type='index').Show()

Set.Plot(Type='adiabatic').Show()