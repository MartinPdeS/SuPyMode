import os
import numpy                 as np
import matplotlib.pyplot     as plt
import matplotlib.gridspec   as gridspec
from numpy                   import min, max
from itertools               import combinations


from SuPyModes.Config        import *
from SuPyModes.Directories   import *
from SuPyModes.utils         import Multipage, prePlot, ToList, Enumerate


class SetPlots(object):

    @prePlot
    def PlotIndex(self, nMax, fig):
        I = self.Index
        for i in range(nMax):
            plt.plot(self.Geometry.ITRList, I[:,i], label=f'{i}')

        return self.PlotKwarg['Index'], [min(I), max(I)]


    @prePlot
    def PlotCoupling(self, nMax, fig):
        C = self.Coupling
        comb = tuple(combinations( np.arange(nMax), 2 ) )

        for n, (i,j) in enumerate( comb ):
            plt.plot(self.Geometry.ITRList[1:], C[:,i,j], label=f'{i} - {j}')

        return self.PlotKwarg['Coupling'], [min(C), max(C)]


    @prePlot
    def PlotAdiabatic(self, nMax, fig):
        A = self.Adiabatic
        comb = tuple(combinations( np.arange(nMax), 2 ) )

        for n, (i,j) in enumerate( comb ):
            plt.plot(self.Geometry.ITRList[1:], A[:,i,j], label=f'{i} - {j}')

        return self.PlotKwarg['Adiabatic'], [min(A), max(A)]


    def PlotFields(self, iter, nMax):

        fig   = plt.figure(figsize=((nMax+1)*3,3))
        spec2 = gridspec.GridSpec(ncols=(nMax+1), nrows=3, figure=fig)

        for mode in range(nMax):
            axes = fig.add_subplot(spec2[0:2,mode])
            str   = r"$n_{eff}$"
            title = f"Mode {mode} [{str}: {self[mode][iter].Index:.6f}]"

            self[mode][iter].__plot__(axes, title)
            fig.suptitle(f'ITR:{self.Geometry.ITRList[iter]}')

        axes = fig.add_subplot(spec2[0:2,-1])

        self.Geometry.__plot__(axes)

        plt.tight_layout()

        return fig


    def GenFigures(self, Input, iter, nMax, PlotKwarg):
        iter = ToList(iter)

        Input = set(Input)

        if not nMax: nMax = len(self.SuperModes)

        self.UpdatePlotKwarg(PlotKwarg)

        figures = []

        if Input & set( ['All', 'Index'] ):     figures.append( self.PlotIndex(nMax=nMax) )

        if Input & set( ['All', 'Coupling'] ):  figures.append( self.PlotCoupling(nMax=nMax) )

        if Input & set( ['All', 'Adiabatic'] ): figures.append( self.PlotAdiabatic(nMax=nMax) )

        for i in iter:
            if Input & set( ['All', 'Fields'] ):    figures.append( self.PlotFields(i, nMax=nMax) )

        return figures


    def Plot(self, Input, iter=0, nMax=None, PlotKwarg=None):
        figures = self.GenFigures(Input, iter, nMax, PlotKwarg)

        plt.show()


    def UpdatePlotKwarg(self, kwargs):
        if isinstance(kwargs, dict):
            for key, val in kwargs.items():
                BasePlotKwarg[key].update(kwargs[key])

        self.PlotKwarg = BasePlotKwarg


    def SaveFig(self, Input, Directory, iter=0, nMax=None, dpi=200, PlotKwarg=None):
        figures = self.GenFigures(Input, iter, nMax, PlotKwarg)

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
            self._Adiabatic = self.CppSolver.ComputingAdiabatic()
            return self._Adiabatic

        else:
            return self._Adiabatic


    @property
    def Betas(self):
        if self._Betas is None:
            self._Betas = self.CppSolver.GetBetas()
            return self._Betas

        else:
            return self._Betas


    @property
    def Index(self):
        if self._Index is None:
            self._Index = self.CppSolver.GetIndices()
            return self._Index

        else:
            return self._Index


    @property
    def M(self):
        if self._M is None:
            self._M = self.Coupling
            for i in range(self._M.shape[2]):
                self._M[i,:,:] += np.diag(self.Beta[i, :] )

            return self._M

        else:
            return self._M















# -
