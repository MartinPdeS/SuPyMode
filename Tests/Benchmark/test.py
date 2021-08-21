from SuPyModes.Geometry          import Geometry, Circle, Fused4, Gradient
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
from SuPyModes.fibers            import *


"""
FIGURE 2.5 SBB_____________________________________________________

"""

A = Fiber_DCF1300S_20(1.55)
B = Fiber_DCF1300S_33(1.55)
C = Fiber_2028M12(1.55)

Clad = Fused3( Radius = 62.5, Fusion = 0.9, Index = Fused_silica(1.55))


Gradient0 = Gradient(Center = Clad.C[0], Nin=A.nClad, Nout = Fused_silica(1.55), Rout=A.rClad*2)
Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad*2, Gradient=Gradient0 )
Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )


Clad1 = Circle( Position = Clad.C[1], Radi = B.rClad, Index = B.nClad )
Core1 = Circle( Position = Clad.C[1], Radi = B.rCore, Index = B.nCore )

Clad2 = Circle( Position = Clad.C[2], Radi = C.rClad, Index = C.nClad )
Core2 = Circle( Position = Clad.C[2], Radi = C.rCore, Index = C.nCore )


SMF28 = Geometry(Objects = [Clad, Clad0, Clad1, Clad2, Core0, Core1, Core2],
                 Xbound  = [-110, 110],
                 Ybound  = [-110, 110],
                 Nx      = 15,
                 Ny      = 15,
                 GConv   = 0)

#SMF28.Plot()

Sol = SuPySolver(Coupler    = SMF28,
                 Tolerance  = 1e-30,
                 MaxIter    = 1000,
                 nMode      = 17,
                 sMode      = 15,
                 Error      = 2)

SuperModes = Sol.GetModes(wavelength    = 1.55,
                          Nstep         = 300,
                          ITRi          = 1,
                          ITRf          = 0.05,
                          RightSymmetry = 0,
                          TopSymmetry   = 0,
                          Sorting       = 'Field')


SuperModes.ComputeM()
