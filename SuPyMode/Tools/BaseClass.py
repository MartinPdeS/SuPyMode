import os
import numpy                 as np
import matplotlib.pyplot     as plt
import matplotlib.gridspec   as gridspec
from numpy                   import min, max
from itertools               import combinations
from cycler                  import cycler
plt.rc('lines', linewidth=2)
plt.rc('axes', prop_cycle=(
                           cycler('linestyle', ['-', '--', ':', '-.']) *
                           cycler('color', ['r', 'g', 'b', 'y', 'k'])
                           ))

from SuPyMode.Tools.Config        import *
from SuPyMode.Tools.Directories   import *
from SuPyMode.Tools.utils         import Multipage, prePlot, ToList, Enumerate




class SetProperties(object):

    @property
    def ITR(self):
        return self.Geometry.ITRList


    @property
    def Coupling(self):
        if self._Coupling is None:
            self._Coupling = self.CppSolver.ComputingCoupling()
            return self._Coupling

        else:
            return self._Coupling


    @property
    def Adiabatic(self):
        if self._Adiabatic is None:
            self._Coupling = self.CppSolver.ComputingCoupling()
            self._Adiabatic = self.CppSolver.ComputingAdiabatic()
            return self._Adiabatic

        else:
            return self._Adiabatic


    @property
    def Beta(self):
        if self._Beta is None:
            self._Beta = self.CppSolver.GetBetas()
            return self._Beta

        else:
            return self._Beta


    @property
    def Index(self):
        if self._Index is None:
            self._Index = self.CppSolver.GetIndices()
            return self._Index

        else:
            return self._Index













# -
