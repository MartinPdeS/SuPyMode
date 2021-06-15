from SuPyModes.utils     import *
from SuPyModes.sellmeier import Fused_silica

class FiberA():
    def __init__(self, wavelength):
        self.nClad = NA2nCore( 0.11, Fused_silica(wavelength)  )
        self.nCore = NA2nCore( 0.12, self.nClad )
        self.rClad = 19.9/2
        self.rCore = 4.6


class FiberB():
    def __init__(self, wavelength):
        self.nClad = NA2nCore( 0.11, Fused_silica(wavelength)  )
        self.nCore = NA2nCore( 0.125, self.nClad )
        self.rClad = 33/2
        self.rCore = 4.5


class SMF28():
    def __init__(self, wavelength):
        self.nClad = Fused_silica(wavelength)
        self.nCore = NA2nCore( 0.14, self.nClad )
        self.rClad = 19.9/2
        self.rCore = 4.6
