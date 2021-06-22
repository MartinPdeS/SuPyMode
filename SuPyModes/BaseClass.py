import os
import logging
import numpy                 as np
import matplotlib.pyplot     as plt
import matplotlib.gridspec   as gridspec
import matplotlib.colors     as colors
from numpy                   import min, max
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools               import combinations
from progressbar             import ProgressBar

from SuPyModes.Config        import *
from SuPyModes.Directories   import *
from SuPyModes.utils         import Multipage, prePlot, ToList, Enumerate


class SetPlots(object):

    @prePlot
    def PlotIndex(self, nMax, fig):
        I = self.Index

        for i in range(nMax):
            plt.plot(self.Geometry.ITRList, I[i], label=f'{i}')

        return self.PlotKwarg['Index'], [min(I), max(I)]


    @prePlot
    def PlotCoupling(self, nMax, fig):
        C = self.Coupling
        comb = tuple(combinations( np.arange(nMax), 2 ) )

        for n, (i,j) in enumerate( comb ):
            plt.plot(self.Geometry.ITRList[1:-1], C[i,j,1:], label=f'{i} - {j}')

        return self.PlotKwarg['Coupling'], [min(C), max(C)]


    @prePlot
    def PlotAdiabatic(self, nMax, fig):
        A = self.Adiabatic
        comb = tuple(combinations( np.arange(nMax), 2 ) )

        for n, (i,j) in enumerate( comb ):
            plt.plot(self.Geometry.ITRList[1:-1], A[i,j,1:], label=f'{i} - {j}')

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
    def Index(self):
        if self._Index is None:
            return self.GetIndex()
        else:
            return self._Index


    @property
    def Beta(self):
        if self._Beta is None:
            return self.GetBeta()
        else:
            return self._Beta


    @property
    def Coupling(self):
        if not isinstance(self._Coupling, np.ndarray):
            return self.GetCoupling()
        else:
            return self._Coupling


    @property
    def Adiabatic(self):
        if not isinstance(self._Adiabatic, np.ndarray):
            return self.GetAdiabatic()
        else:
            return self._Adiabatic


    @property
    def M(self):
        if self._M is None:
            self._M = self.Coupling[:,:]
            for i in range(self._M.shape[2]):
                self._M[:,:,i] += np.diag(self.Beta[:, i] )

            return self._M

        else:
            return self._M


    def GetIndex(self):
        logging.info('Computing effective indices...')
        Index = []

        for i, Supermode in enumerate( self.SuperModes ):
            for j, iter in enumerate( Supermode.Slice ):
                Index.append(iter.Index)

        Index = np.asarray( Index ).reshape([i+1, j+1])

        return Index


    def GetBeta(self):
        logging.info('Computing propagation constant...')
        Beta = []

        for i, Supermode in enumerate( self.SuperModes ):
            for j, iter in enumerate( Supermode.Slice ):
                Beta.append(iter.Beta)

        Beta = np.asarray( Beta ).reshape([i+1, j+1])

        return Beta


    def GetCoupling(self):
        C = []
        for n, (i,j) in Enumerate( self.combinations, msg='Computing mode coupling... ' ):

            C.append( self[i].GetCoupling(self[j]) )


        tri = np.zeros( ( self.NSolutions, self.NSolutions, len(self.Geometry.ITRList)-1 ) )
        tri[np.triu_indices(self.NSolutions, 1)] = C
        tri[np.tril_indices(self.NSolutions, -1)] = C


        self._Coupling = tri

        return self._Coupling


    def GetAdiabatic(self):
        A = []
        for n, (i,j) in Enumerate( self.combinations, msg='Computing adiabatic criterion... '):

            A.append( self[i].GetAdiabatic(self[j]) )

        tri = np.zeros( ( self.NSolutions, self.NSolutions, len(self.Geometry.ITRList)-1 ) )
        tri[np.triu_indices(self.NSolutions, 1)]  = A
        tri[np.tril_indices(self.NSolutions, -1)] = A

        self._Adiabatic = tri

        return self._Adiabatic
















# -
