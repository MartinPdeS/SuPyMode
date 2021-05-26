import os
import logging
import numpy                 as np
import matplotlib.pyplot     as plt
import matplotlib.gridspec   as gridspec
import matplotlib.colors     as colors
from numpy                   import min, max
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools               import combinations

from SuPyModes.Config        import *
from SuPyModes.utils         import RecomposeSymmetries, Multipage, prePlot
from SuPyModes.Directories   import *


class SetPlots(object):

    @prePlot
    def PlotIndex(self, nMax, fig):
        I = self.Index

        for i in range(nMax):
            plt.plot(self.ITR, I[i], label=f'{i}')

        return self.PlotKwarg['Index'], [min(I), max(I)]


    @prePlot
    def PlotCoupling(self, nMax, fig):
        C = self.Coupling
        comb = tuple(combinations( np.arange(nMax), 2 ) )

        for n, (i,j) in enumerate( comb ):
            plt.plot(self.ITR[1:-1], C[i,j,1:], label=f'{i} - {j}')

        return self.PlotKwarg['Coupling'], [min(C), max(C)]


    @prePlot
    def PlotAdiabatic(self, nMax, fig):
        A = self.Adiabatic
        comb = tuple(combinations( np.arange(nMax), 2 ) )

        for n, (i,j) in enumerate( comb ):
            plt.plot(self.ITR[1:-1], A[i,j,1:], label=f'{i} - {j}')

        return self.PlotKwarg['Adiabatic'], [min(A), max(A)]


    def PlotFields(self, iter, nMax):

        fig   = plt.figure(figsize=((nMax+1)*3,3))
        spec2 = gridspec.GridSpec(ncols=(nMax+1), nrows=3, figure=fig)

        for mode in range(nMax):
            Field, xaxis, yaxis = self[mode].FullField(iter)
            axes = fig.add_subplot(spec2[0:2,mode])
            axes.pcolormesh(yaxis, xaxis, Field, shading='auto')
            axes.set_aspect('equal')
            str   = r"$n_{eff}$"
            title = f"Mode {mode} [{str}: {self[mode].Index[iter]:.6f}]"
            axes.set_title(title, fontsize=8)
            if mode == 0:
                axes.set_ylabel(r'Y-distance [$\mu$m]', fontsize=6)
                axes.set_xlabel(r'X-distance [$\mu$m]', fontsize=6)

        axes = fig.add_subplot(spec2[0:2,-1])

        Profile, _, _ = RecomposeSymmetries(Input      = self.IndexProfile.T,
                                            Symmetries = self.Symmetries,
                                            Xaxis      = None,
                                            Yaxis      = None)

        axes.set_title('Rasterized RI profil', fontsize=10)
        axes.set_ylabel(r'Y-distance [$\mu$m]')
        axes.set_xlabel(r'X-distance [$\mu$m]')

        Indices = np.unique(np.abs(Profile))
        vmin=sorted(Indices)[1]/1.1
        vmax=sorted(Indices)[-1]

        pcm = axes.pcolormesh( yaxis,
                               xaxis,
                               np.abs(Profile),
                               cmap    = plt.cm.coolwarm,
                               norm=colors.LogNorm(vmin=vmin, vmax=vmax),
                               shading='auto'
                              )

        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=colors.LogNorm(vmin=vmin, vmax=vmax))

        axes.contour(yaxis, xaxis, np.abs(Profile), levels=Indices, colors='k')

        sm._A = []
        cbar = plt.colorbar(sm, ax=axes, cax=cax)

        axes.set_aspect('equal')

        plt.tight_layout()

        return fig



    def GenFigures(self, Input, iter, nMax, PlotKwarg):
        Input = set(Input)

        if not nMax: nMax = len(self.SuperModes)

        self.UpdatePlotKwarg(PlotKwarg)

        figures = []

        if Input & set( ['All', 'Index'] ):     figures.append( self.PlotIndex(nMax=nMax) )

        if Input & set( ['All', 'Coupling'] ):  figures.append( self.PlotCoupling(nMax=nMax) )

        if Input & set( ['All', 'Adiabatic'] ): figures.append( self.PlotAdiabatic(nMax=nMax) )

        if Input & set( ['All', 'Fields'] ):    figures.append( self.PlotFields(iter, nMax=nMax) )

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
    def Index(self):
        logging.info('Computing effective indices...')
        I = []
        for i in range(self.NSolutions):
            I.append( self[i].Index )

        return I


    @property
    def Beta(self):
        logging.info(r'Computing mode propagation constant ($\beta$)...')
        B = []
        for i in range(self.NSolutions):
            B.append( self[i].Beta )

        return B


    @property
    def Coupling(self):
        if not isinstance(self._Coupling, np.ndarray):
            self.GetCoupling()
            return self._Coupling
        else:
            return self._Coupling


    @property
    def Adiabatic(self):
        if not isinstance(self._Adiabatic, np.ndarray):
            self.GetAdiabatic()
            return self._Adiabatic
        else:
            return self._Adiabatic


    def GetCoupling(self):
        logging.info('Computing mode coupling...')
        C = []
        for (i,j) in self.combinations:
            C.append( self[i].GetCoupling(self[j]) )

        tri = np.zeros( ( self.NSolutions, self.NSolutions, len(self.ITR)-1 ) )
        tri[np.triu_indices(self.NSolutions, 1)] = C
        tri[np.tril_indices(self.NSolutions, -1)] = C


        self._Coupling = tri

        return tri


    def GetAdiabatic(self):
        logging.info('Computing adiabatic criterion...')
        A = []
        for (i,j) in self.combinations:
            A.append( self[i].GetAdiabatic(self[j]) )

        tri = np.zeros( ( self.NSolutions, self.NSolutions, len(self.ITR)-1 ) )
        tri[np.triu_indices(self.NSolutions, 1)]  = A
        tri[np.tril_indices(self.NSolutions, -1)] = A

        self._Adiabatic = tri

        return tri
