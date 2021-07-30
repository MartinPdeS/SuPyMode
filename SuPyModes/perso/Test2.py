from SuPyModes.Geometry          import Geometry, Circle, Fused2
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica

from SuPyModes.perso.fibers import *
import numpy as np
import matplotlib.pyplot as plt




"""
FIGURE 2.5 SBB_____________________________________________________

"""
wavelength = 1.55
k          = 2*np.pi/wavelength
Nstep      = 300
sMode      = 12
IndexSilica = Fused_silica(wavelength)
Clad = Circle( Position = (0,0), Radi = 62.5, Index =  IndexSilica)
Core0 = Circle( Position = (0,0), Radi = 4.1, Index = IndexSilica+0.005 )



SMF_silica = Geometry(Objects = [Clad, Core0],
                       Xbound  = [-90, 90],
                       Ybound  = [-90, 90],
                       Nx      = 120,
                       Ny      = 120)

SMF_silica.Plot()

Sol = SuPySolver(Coupler    = SMF_silica,
                 Tolerance  = 1e-30,
                 MaxIter    = 1000,
                 nMode      = sMode+2,
                 sMode      = sMode,
                 Error      = 2)

SuperModes = Sol.GetModes(wavelength    = wavelength,
                          Nstep         = Nstep,
                          ITRi          = 1,
                          ITRf          = 0.05,
                          RightSymmetry = 0,
                          TopSymmetry   = 0,
                          Sorting       = 'Field')

SuperModes.Plot(['Fields', 'Adiabatic', 'Index'])
