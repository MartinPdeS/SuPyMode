import os
import logging
import numpy               as np
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
from itertools             import combinations

from SuPyModes.Config      import *
from SuPyModes.utils       import RecomposeSymmetries, Multipage
from SuPyModes.Directories import *

class SetPlots(object):
    def GetFigure(self, Dict, kwargs):
        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(111)
        ax.grid()
        ax.set_ylabel(Dict['name'] + Dict['unit'])
        ax.set_xlabel('ITR')
        if kwargs['xlim'] is not None: ax.set_xlim(left  = kwargs['xlim'][0],
                                                   right = kwargs['xlim'][1])

        if kwargs['ylim'] is not None: ax.set_ylim( bottom = kwargs['ylim'][0],
                                                    top    = kwargs['ylim'][1] )

        if kwargs['yscale'] == 'log': ax.set_yscale('log')
        if kwargs['xscale'] == 'log': ax.set_xscale('log')
        return fig


    def PlotIndex(self, nMax):
        I = self.Index
        fig = self.GetFigure(IndexDict, self.PlotKwarg['Index'])

        for i in range(nMax):
            plt.plot(self.ITR, I[i], label=f'{i}')
        plt.legend(fontsize=6, title="Modes", fancybox=True)

        return fig


    def PlotCoupling(self, nMax):
        C = self.Coupling
        comb = tuple(combinations( np.arange(nMax), 2 ) )
        fig = self.GetFigure(CouplingDict, self.PlotKwarg['Coupling'])

        for n, (i,j) in enumerate( comb ):
            plt.plot(self.ITR[1:-1], C[i,j,1:], label=f'{i} - {j}')
        plt.legend(fontsize=6, title="Modes", fancybox=True)

        return fig


    def PlotAdiabatic(self, nMax):
        A = self.Adiabatic
        fig = self.GetFigure(AdiabaticDict, self.PlotKwarg['Adiabatic'])
        comb = tuple(combinations( np.arange(nMax), 2 ) )

        for n, (i,j) in enumerate( comb ):
            plt.semilogy(self.ITR[1:-1], A[i,j,1:], label=f'{i} - {j}')
        plt.legend(fontsize=6, title="Modes", fancybox=True)

        return fig


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

        axes.pcolormesh(yaxis, xaxis, np.abs(Profile), shading='auto')
        axes.set_aspect('equal')
        title = f"Index profile [ITR: {self.ITR[iter]:.3f}]"
        axes.set_title(title, fontsize=8)
        plt.tight_layout()

        return fig



    def Plot(self, Input, iter=0, nMax=None, PlotKwarg=None):
        if not nMax: nMax = len(self.SuperModes)

        self.UpdatePlotWkargs(PlotKwargs)

        figures = []

        if 'Index'     in Input: figures.append( self.PlotIndex(nMax=nMax, kwargs=kwargs) )

        if 'Coupling'  in Input: figures.append( self.PlotCoupling(nMax=nMax, kwargs=kwargs) )

        if 'Adiabatic' in Input: figures.append( self.PlotAdiabatic(nMax=nMax, kwargs=kwargs) )

        if 'Fields'    in Input: figures.append( self.PlotFields(iter, nMax=nMax, kwargs=kwargs) )

        plt.show()


    def UpdatePlotKwarg(self, kwargs):
        if isinstance(kwargs, dict):
            for key, val in kwargs.items():
                BasePlotKwarg[key].update(kwargs[key])

        self.PlotKwarg = BasePlotKwarg


    def Save(self, Input, Directory, iter=0, nMax=None, dpi=200, PlotKwarg=None):

        if not nMax: nMax = len(self.SuperModes)

        self.UpdatePlotKwarg(PlotKwarg)

        figures = []

        if 'Index'     in Input: figures.append( self.PlotIndex(nMax=nMax) )

        if 'Coupling'  in Input: figures.append( self.PlotCoupling(nMax=nMax) )

        if 'Adiabatic' in Input: figures.append( self.PlotAdiabatic(nMax=nMax) )

        if 'Fields'    in Input: figures.append( self.PlotFields(iter, nMax=nMax) )

        dir = os.path.join(RootPath, Directory)

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
