from SuPyModes.utils     import *
from SuPyModes.sellmeier import Fused_silica
from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica


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


class Fiber_SMF28():
    def __init__(self, wavelength):
        self.nClad = Fused_silica(wavelength)
        self.nCore = NA2nCore( 0.14, self.nClad )
        self.rClad = 19.9/2
        self.rCore = 4.1





A, B = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = Fused_silica(1.55)  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


BBA_silica = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-150, 150],
               Ybound  = [-150, 150],
               Nx      = 150,
               Ny      = 150)

BBA_silica.Rotate(28)


B, A = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = Fused_silica(1.55)  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


AAB_silica = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-150, 150],
               Ybound  = [-150, 150],
               Nx      = 150,
               Ny      = 150)

AAB_silica.Rotate(28)




A, B = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = 1.433  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


BBA_fluoride = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
               Xbound  = [-150, 150],
               Ybound  = [-150, 150],
               Nx      = 150,
               Ny      = 150)

BBA_fluoride.Rotate(28)


B, A = Fiber_DCF1300S_33(wavelength=1.55), Fiber_DCF1300S_20(wavelength=1.55)


Capillary = Circle( Position = [0,0], Radi = 140, Index = 1.433  )

Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))

Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad, Index = A.nClad )
Clad1 = Circle( Position = Clad.C[1], Radi = A.rClad, Index = A.nClad )
Clad2 = Circle( Position = Clad.C[2], Radi = B.rClad, Index = B.nClad )

Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )
Core1 = Circle( Position = Clad.C[1], Radi = A.rCore, Index = A.nCore )
Core2 = Circle( Position = Clad.C[2], Radi = B.rCore, Index = B.nCore )


AAB_fluoride = Geometry(Objects = [Capillary, Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
                       Xbound  = [-150, 150],
                       Ybound  = [-150, 150],
                       Nx      = 150,
                       Ny      = 150)

AAB_fluoride.Rotate(28)
