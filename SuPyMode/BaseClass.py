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

from SuPyMode.Config        import *
from SuPyMode.Directories   import *
from SuPyMode.utils         import Multipage, prePlot, ToList, Enumerate


class SetPlots(object):

    @prePlot
    def PlotIndex(self, fig):
        I = self.Index
        for i in range(self.sMode):
            plt.plot(self.Geometry.ITRList, I[:,i], label=f'{i}')

        return self.PlotKwarg['Index'], [min(I), max(I)]

    @prePlot
    def PlotBeta(self, fig):
        B = self.Beta
        for i in range(self.sMode):
            plt.plot(self.Geometry.ITRList, B[:,i], label=f'{i}')

        return self.PlotKwarg['Beta'], [min(B), max(B)]


    @prePlot
    def PlotCoupling(self, fig):
        C = self.Coupling

        for n, (i,j) in enumerate( self.Combination ):
            plt.plot(self.Geometry.ITRList[0:], C[:,i,j], label=f'{i} - {j}')

        return self.PlotKwarg['Coupling'], [min(C), max(C)]


    @prePlot
    def PlotAdiabatic(self, fig):
        A = self.Adiabatic
        print(A)
        print('='*100)
        print(self.Coupling)
        for n, (i,j) in enumerate( self.Combination ):
            plt.plot(self.Geometry.ITRList[0:], A[:,i,j], label=f'{i} - {j}')

        return self.PlotKwarg['Adiabatic'], [min(A), max(A)]


    def PlotFields(self, iter):

        fig   = plt.figure(figsize=((self.sMode+1)*3,3))
        spec2 = gridspec.GridSpec(ncols=(self.sMode+1), nrows=3, figure=fig)

        for mode in range(self.sMode):
            axes = fig.add_subplot(spec2[0:2,mode])
            str   = r"$n_{eff}$"
            title = f"Mode {mode} [{str}: {self[mode][iter].Index:.6f}]"

            self[mode][iter].__plot__(axes, title)
            fig.suptitle(f'ITR:{self.Geometry.ITRList[iter]}')

        axes = fig.add_subplot(spec2[0:2,-1])

        self.Geometry.__plot__(axes)

        plt.tight_layout()

        return fig


    def GenFigures(self, Input, iter, PlotKwarg):
        iter = ToList(iter)

        Input = set(Input)

        self.UpdatePlotKwarg(PlotKwarg)

        figures = []

        if Input & set( ['All', 'Index'] ):     figures.append( self.PlotIndex() )

        if Input & set( ['All', 'Beta'] ):     figures.append( self.PlotBeta() )

        if Input & set( ['All', 'Coupling'] ):  figures.append( self.PlotCoupling() )

        if Input & set( ['All', 'Adiabatic'] ): figures.append( self.PlotAdiabatic() )

        for i in iter:
            if Input & set( ['All', 'Fields'] ):    figures.append( self.PlotFields(i) )

        return figures


    def Plot(self, Input, iter=0, PlotKwarg=None, Combination=None):

        if Combination is None:
            self.Combination = tuple(combinations( np.arange(self.sMode), 2 ) )
        else:
            self.Combination = Combination

        figures = self.GenFigures(Input, iter, PlotKwarg)

        plt.show()


    def UpdatePlotKwarg(self, kwargs):
        if isinstance(kwargs, dict):
            for key, val in kwargs.items():
                BasePlotKwarg[key].update(kwargs[key])

        self.PlotKwarg = BasePlotKwarg


    def SaveFig(self, Input, Directory, iter=0, dpi=100, PlotKwarg=None, Combination=None):

        if Combination is None:
            self.Combination = tuple(combinations( np.arange(self.sMode), 2 ) )
        else:
            self.Combination = Combination

        figures = self.GenFigures(Input, iter, PlotKwarg)

        dir = os.path.join(ZeroPath, Directory) + '.pdf'

        Multipage(dir, figs=figures, dpi=dpi)



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
