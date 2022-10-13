"""
3x3 Coupler
===========

.. image:: ../../images/Example3/Geometry.png
   :width: 400
   :align: center

.. image:: ../../images/Example3/Fields.png
   :width: 800
   :align: center


.. image:: ../../images/Example3/Index.png
   :width: 600
   :align: center


.. image:: ../../images/Example3/Adiabatic.png
   :width: 600
   :align: center
"""

# sphinx_gallery_thumbnail_path = '../images/Example3/Geometry.png'



from FiberFusing                import Geometry, Fused3, Circle, BackGround
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica

Air = BackGround(Index=1) 

Wavelength = 1.55e-6
Index = FusedSilica.GetRI(Wavelength)

Clad = Fused3(FiberRadius = 60, Fusion = 0.3, Index = Index)

Cores = [ Circle(Position=Core, Radius=4.1, Index=Index+0.005) for Core in Clad.Cores]

Geo = Geometry(Objects = [Air, Clad] + Cores,
               Xbound  = [-120, 120],
               Ybound  = [-100, 130],
               Nx      = 180,
               Ny      = 180)

Geo.Plot().Show()

Sol = SuPySolver(Geometry=Geo, Tolerance=1e-8, MaxIter = 10000)

Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)

Sol.AddModes(Sorting         = 'Index',
             Symmetries      = {'Right': 0, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nComputedMode   = 6,
             nSortedMode     = 4 )

Set = Sol.GetSet()

Set.Plot(Type='field', ITR=[1, 0.3]).Show()

Set.Plot(Type='index').Show()

Set.Plot(Type='adiabatic').Show()