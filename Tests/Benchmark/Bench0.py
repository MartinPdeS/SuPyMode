from SuPyModes.Geometry          import Geometry, Circle, Fused2
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica



"""
FIGURE 2.5 SBB_____________________________________________________

"""

Clad = Circle( Position = (0,0),
               Radi     = 62.5,
               Index    =  Fused_silica(1.55))

Core0 = Circle( Position = (0,0),
                Radi     = 4.1,
                Index    = Fused_silica(1.55)+0.005 )

SMF28 = Geometry(Objects = [Clad, Core0],
                 Xbound  = [-90, 90],
                 Ybound  = [-90, 90],
                 Nx      = 120,
                 Ny      = 120,
                 GConv   = 1.0)

SMF28.Plot()

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

SuperModes.Plot( Input=['Fields', 'Adiabatic'],
                 Combination = [(0,1), (0,10)] )





"""
Note that the bad correspondence are due to the fact that we do not break the symmetries and
as such the modes are rotating from one slice to the other!

"""
