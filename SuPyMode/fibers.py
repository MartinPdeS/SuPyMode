from SuPyMode.utils     import *
from SuPyMode.sellmeier import Fused_silica
from SuPyMode.Geometry          import Geometry, Circle, Fused3
from SuPyMode.Solver            import SuPySolver
from SuPyMode.sellmeier         import Fused_silica


class Fiber_DCF1300S_20():
    def __init__(self, wavelength):
        self.nClad = NA2nCore( 0.11, Fused_silica(wavelength)  )
        self.nCore = NA2nCore( 0.12, self.nClad )
        self.rClad = 19.9/2
        self.rCore = 4.6


class Fiber_DCF1300S_33():
    def __init__(self, wavelength):
        self.nClad = NA2nCore( 0.11, Fused_silica(wavelength)  )
        self.nCore = NA2nCore( 0.125, self.nClad )
        self.rClad = 33/2
        self.rCore = 4.5


class Fiber_2028M24():
    def __init__(self, wavelength):
        self.nClad = NA2nCore( 0.19, Fused_silica(wavelength)  )
        self.nCore = NA2nCore( 0.11, self.nClad )
        self.rClad = 14.1/2
        self.rCore = 2.3/2

class Fiber_2028M21():
    def __init__(self, wavelength):
        self.nClad = NA2nCore( 0.19, Fused_silica(wavelength)  )
        self.nCore = NA2nCore( 0.11, self.nClad )
        self.rClad = 17.6/2
        self.rCore = 2.8/2


class Fiber_2028M12():
    def __init__(self, wavelength):
        self.nClad = NA2nCore( 0.19, Fused_silica(wavelength)  )
        self.nCore = NA2nCore( 0.11, self.nClad )
        self.rClad = 25.8/2
        self.rCore = 4.1/2


class Fiber_SMF28():
    def __init__(self, wavelength):
        self.nClad = Fused_silica(wavelength)
        self.nCore = Fused_silica(wavelength)+0.005#NA2nCore( 0.14, self.nClad )
        self.rClad = 19.9/2
        self.rCore = 4.1
