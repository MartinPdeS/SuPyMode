import os
import numpy                 as np
import matplotlib.pyplot     as plt

from itertools               import combinations
from cycler                  import cycler
plt.rc('lines', linewidth=2)
plt.rc('axes', prop_cycle=(
                           cycler('linestyle', ['-', '--', ':', '-.']) *
                           cycler('color', ['r', 'g', 'b', 'y', 'k'])
                           ))


from SuPyMode.Tools.Directories   import ZeroPath
from SuPyMode.Tools.utils         import Multipage, ToList
from SuPyMode.Plotting.PlotsUtils import FieldMap
from SuPyMode.Plotting.Plots      import Scene2D
from SuPyMode.Plotting.Plots      import Scene2D, Scene, Axis, Line, Mesh


class SetProperties(object):
    _Coupling     = None
    _Adiabatic    = None
    _Index        = None
    _Beta         = None
    _M            = None
    SuperModes    = []

    @property
    def Geometry(self):
        return self.ParentSolver.Geometry

    @property
    def ITR(self):
        return self.Geometry.ITRList

    @property
    def Symmetries(self):
        return self.Geometry.Symmetries




class SetPlottings():


    def GetCombinations(self):
        Combination = []
        for Mode0 in self.SuperModes:
            for Mode1 in self.SuperModes:
                if Mode0.SolverNumber != Mode1.SolverNumber: continue

                if Mode0.BindingNumber == Mode1.BindingNumber: continue

                if (Mode1, Mode0) in Combination: continue

                Combination.append((Mode0, Mode1))

        return Combination


    def PlotIndex(self, ReturnFig=False):
        Fig = Scene('SuPyMode Figure', UnitSize=(10,4))

        ax = Axis(Row    = 0,
                  Col    = 0,
                  xLabel = 'ITR',
                  yLabel = r'Effective refraction index',
                  Title  = None,
                  Grid   = True,
                  xScale = 'linear',
                  yScale = 'linear')

        for supermode in self.SuperModes:
            supermode._PlotIndex(ax)

        Fig.AddAxes(ax)

        if ReturnFig:
            Fig.Render()
            return Fig.Figure

        Fig.Show()


    def PlotBetas(self, ReturnFig=False):
        Fig = Scene('SuPyMode Figure', UnitSize=(10,4))

        ax = Axis(Row    = 0,
                  Col    = 0,
                  xLabel = 'ITR',
                  yLabel = r'Propagation constante $\beta$',
                  Title  = None,
                  Grid   = True,
                  xScale = 'linear',
                  yScale = 'linear')

        for supermode in self.SuperModes:
            supermode._PlotBetas(ax)

        Fig.AddAxes(ax)

        if ReturnFig:
            Fig.Render()
            return Fig.Figure

        Fig.Show()



    def PlotCoupling(self, ReturnFig=False):
        Fig = Scene('SuPyMode Figure', UnitSize=(10,4))

        Combination = self.GetCombinations()

        for Mode0, Mode1 in Combination:
                ax = Axis(Row    = 0,
                          Col    = 0,
                          xLabel = 'ITR',
                          yLabel = r'Effective refraction index',
                          Title  = None,
                          Grid   = True,
                          xScale = 'linear',
                          yScale = 'linear')

                artist = Line(X=self.ITRList, Y=Mode0.GetCoupling(Mode1), Label=f'{Mode0.Name} - {Mode1.Name}')

                ax.AddArtist(artist)

                Fig.AddAxes(ax)

        if ReturnFig:
            Fig.Render()
            return Fig.Figure

        Fig.Show()


    def PlotAdiabatic(self, ReturnFig=False):
        Fig = Scene('SuPyMode Figure', UnitSize=(10,4))

        Combination = self.GetCombinations()

        for Mode0, Mode1 in Combination:
                ax = Axis(Row     = 0,
                          Col     = 0,
                          xLabel  = 'ITR',
                          yLabel  = r'Adiabatic criterion',
                          Title   = None,
                          Grid    = True,
                          xScale  = 'linear',
                          yLimits = [1e-6, 1e-1],
                          yScale  = 'log')

                artist = Line(X=self.ITRList, Y=Mode0.Binded.GetAdiabaticSpecific(Mode1.Binded), Label=f'{Mode0.Name} - {Mode1.Name}')

                ax.AddArtist(artist)

                Fig.AddAxes(ax)

        if ReturnFig:
            Fig.Render()
            return Fig.Figure

        Fig.Show()


    def PlotFields(self, Slices=[0], ReturnFig=False):

        Fig = Scene('SuPyMode Figure', UnitSize=(3,3))

        for m, Mode in enumerate(self.SuperModes):
            for n, slice in enumerate(Slices):
                ax = Axis(Row              = n,
                          Col              = m,
                          xLabel           = r'',
                          yLabel           = r'',
                          Title            = f'{self[m].Name}  [ITR: {self.ITRList[slice]:.2f}]',
                          Legend           = False,
                          ColorBar         = True,
                          ColorbarPosition = 'right',
                          Grid             = True,
                          Equal            = True,
                          xScale           = 'linear',
                          yScale           = 'linear')

                self[m]._PlotFields(ax, slice)

                Fig.AddAxes(ax)


        if ReturnFig:
            Fig.Render()
            return Fig.Figure

        Fig.Show()





    def ExportPDF(self, Directory='ExportedPDF', Slices=[0, -1], dpi=100):
        Fig0 = self.PlotIndex(ReturnFig=True)
        Fig1 = self.PlotAdiabatic(ReturnFig=True)
        Fig2 = self.PlotFields(Slices=Slices, ReturnFig=True)

        dir = os.path.join(ZeroPath, Directory) + '.pdf'

        Multipage(dir, figs=[Fig0, Fig1, Fig2], dpi=dpi)











# -
