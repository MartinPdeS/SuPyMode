from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

from SuPyModes.perso.fibers import Fiber_DCF1300S_33, Fiber_DCF1300S_20, Fiber_SMF28, AAB_fluoride, BBA_fluoride, AAB_silica, BBA_silica

AAB_fluoride.Plot()
BBA_fluoride.Plot()

Sol = SuPySolver(Coupler    = AAB_fluoride,
                 Tolerance  = 1e-20,
                 MaxIter    = 1000,
                 nMode      = 7,
                 sMode      = 5)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 300,
                          ITRi       = 1,
                          ITRf       = 0.05)


SuperModes.SaveFig(Directory  = 'BBA_fluoride',
                   Input      = ['All'],
                   nMax       = 5)





Sol = SuPySolver(Coupler    = BBA_fluoride,
                 Tolerance  = 1e-20,
                 MaxIter    = 1000,
                 nMode      = 7,
                 sMode      = 5)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 300,
                          ITRi       = 1,
                          ITRf       = 0.05)


SuperModes.SaveFig(Directory  = 'AAB_fluoride',
                   Input      = ['All'],
                   nMax       = 5)







Sol = SuPySolver(Coupler    = AAB_silica,
                 Tolerance  = 1e-20,
                 MaxIter    = 1000,
                 nMode      = 7,
                 sMode      = 5)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 300,
                          ITRi       = 1,
                          ITRf       = 0.05)


SuperModes.SaveFig(Directory  = 'BBA_silica',
                   Input      = ['All'],
                   nMax       = 5)





Sol = SuPySolver(Coupler    = BBA_silica,
                 Tolerance  = 1e-20,
                 MaxIter    = 1000,
                 nMode      = 7,
                 sMode      = 5)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 300,
                          ITRi       = 1,
                          ITRf       = 0.05)


SuperModes.SaveFig(Directory  = 'AAB_silica',
                   Input      = ['All'],
                   nMax       = 5)
