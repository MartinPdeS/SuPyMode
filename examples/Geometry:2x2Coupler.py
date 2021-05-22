from SuPyModes.Geometry          import Geometry, Fused2
from SuPyModes.sellmeier         import Fused_silica

Clad = Fused2(Radius =  62.5,
              Fusion  = 1,
              Index   = Fused_silica(1.55))

Clad.Plot()
