from SuPyModes.Geometry          import Geometry, Circle, Fused2
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica



"""
FIGURE 2.5 SBB_____________________________________________________

"""

Clad = Fused2(Radius  =  62.5,
              Fusion  = 0.9,
              Index   = Fused_silica(1.55))

Core0 = Circle( Position = Clad.C[0],
                Radi     = 4.1,
                Index    = Fused_silica(1.55)+0.005 )

Core1 = Circle( Position = Clad.C[1],
                Radi     = 4.1,
                Index    = Fused_silica(1.55)+0.005 )

Coupler = Geometry(Objects = [Clad, Core0, Core1],
                   Xbound  = [-120, 0],
                   Ybound  = [-120, 0],
                   Nx      = 120,
                   Ny      = 120,
                   GConv   = 0.0)

Coupler.Plot()

Sol = SuPySolver(Coupler    = Coupler,
                 Tolerance  = 1e-30,
                 MaxIter    = 1000,
                 nMode      = 6,
                 sMode      = 5,
                 Error      = 2)



SuperModes = Sol.GetModes(wavelength    = 1.55,
                          Nstep         = 500,
                          ITRi          = 1,
                          ITRf          = 0.05,
                          RightSymmetry = 1,
                          TopSymmetry   = 1,
                          Sorting       = 'Field')

SuperModes.Plot( Input=['Fields', 'Adiabatic'],
                 iter = [0,-1],
                 Combination = [(0,1), (0,2), (0,3)]
                 )


SuperModes = Sol.GetModes(wavelength    = 1.55,
                          Nstep         = 500,
                          ITRi          = 1,
                          ITRf          = 0.05,
                          RightSymmetry = -1,
                          TopSymmetry   = 1,
                          Sorting       = 'Field')

SuperModes.Plot( Input=['Fields', 'Adiabatic'],
                 iter = [0,-1],
                 Combination = [(0,1), (0,2), (0,3)]
                 )
