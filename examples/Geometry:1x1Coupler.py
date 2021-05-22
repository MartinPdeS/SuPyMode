from SuPyModes.Geometry          import Geometry, Circle
from SuPyModes.sellmeier         import Fused_silica

Clad = Circle(Radi = 62.5, Position = (0,0), Index = Fused_silica(1.55))

Clad.Plot()
