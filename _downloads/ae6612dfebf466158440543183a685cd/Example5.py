"""
7x7 Coupler
===========

.. image:: ../../images/Example5/Geometry.png
   :width: 400
   :align: center

.. image:: ../../images/Example5/Fields.png
   :width: 800
   :align: center


.. image:: ../../images/Example5/Index.png
   :width: 600
   :align: center


.. image:: ../../images/Example5/Adiabatic.png
   :width: 600
   :align: center
"""

# sphinx_gallery_thumbnail_path = '../images/Example5/Geometry.png'



from FiberFusing                import Geometry, Fused7, Circle, BackGround
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica

Air = BackGround(Index=1) 

Wavelength = 1.55e-6
Index = FusedSilica.GetRI(Wavelength)

Clad = Fused7(FiberRadius=62.5, Fusion=0.3, Index=Index)

Cores = [ Circle(Position=Core, Radius=4.1, Index=Index+0.005) for Core in Clad.Cores]

Geo = Geometry(Objects = [Air, Clad] + Cores,
               Xbound  = [-190, 0],
               Ybound  = [-190, 190],
               Nx      = 50,
               Ny      = 100)

# Geo.Plot().Show()

Sol = SuPySolver(Geometry=Geo, Tolerance=1e-8, MaxIter=10000)

Sol.CreateSuperSet(Wavelength=1.55, NStep=10, ITRi=1, ITRf=0.99)


Sol.AddModes(Sorting         = 'Field',
             Symmetries      = {'Right': 1, 'Left': 0, 'Top': 0, 'Bottom': 0},
             nComputedMode   = 30,
             nSortedMode     = 2 )

# Sol.AddModes(Sorting         = 'Field',
#              Symmetries      = {'Right': -1, 'Left': 0, 'Top': 0, 'Bottom': 0},
#              nComputedMode   = 6,
#              nSortedMode     = 2 )

Set = Sol.GetSet()


Set.Plot(Type='field', Slice=[1, 3, 4, 5]).Show()

Set.Plot(Type='index').Show()

Set.Plot(Type='adiabatic').Show()