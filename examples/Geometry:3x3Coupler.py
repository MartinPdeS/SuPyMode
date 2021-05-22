from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.sellmeier         import Fused_silica

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad.Plot()
