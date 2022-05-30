
from SuPyMode.Tools.utils       import NA2nCore
from SuPyMode.Geometry          import Geometry, Circle, Fused3
from SuPyMode.Solver            import SuPySolver
from PyOptik                    import ExpData
from SuPyMode.Geometry          import Circle

class Fiber_DCF1300S_20():
    def __init__(self, Wavelength):
        Index = ExpData('FusedSilica').GetRI(Wavelength*1e-6)
        self.nClad = NA2nCore( 0.11, Index )
        self.nCore = NA2nCore( 0.12, self.nClad )
        self.rClad = 19.9/2
        self.rCore = 4.6

    def Get(self, Position):
        self.Fiber = [
                       Circle( Position=Position, Radius=self.rClad, Index=self.nClad ),
                       Circle( Position=Position, Radius=self.rCore, Index=self.nCore ),
                       ]
        return self.Fiber


class Fiber_DCF1300S_33():
    def __init__(self, Wavelength):
        Index = ExpData('FusedSilica').GetRI(Wavelength*1e-6)
        self.nClad = NA2nCore( 0.11, Index  )
        self.nCore = NA2nCore( 0.125, self.nClad )
        self.rClad = 33/2
        self.rCore = 4.5

    def Get(self, Position):
        self.Fiber = [
                       Circle( Position=Position, Radius=self.rClad, Index=self.nClad ),
                       Circle( Position=Position, Radius=self.rCore, Index=self.nCore ),
                       ]
        return self.Fiber


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
