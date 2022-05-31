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



class SetPlottings():

    def PlotIndex(self, Scene, Col, Row):
        for supermode in self.SuperModes:
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.ITRList,
                          y        = supermode.Index,
                          Fill     = False,
                          Legend   = supermode.ModeNumber,
                          xLabel   = r'ITR',
                          yLabel   = r'Propagation constante $\beta$')

            Scene.SetAxes(Col, Row, Equal=False, Legend=True, yLimits=[self.Geometry.MinIndex/1.005, self.Geometry.MaxIndex], LegendTitle='Mode')


    def PlotBeta(self, Scene, Col, Row):
        for supermode in self.SuperModes:
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.ITRList,
                          y        = supermode.Beta,
                          Fill     = False,
                          Legend   = supermode.ModeNumber,
                          xLabel   = r'ITR',
                          yLabel   = r'Propagation constante $\beta$')

            Scene.SetAxes(Col, Row, Equal=False, Legend=True, LegendTitle='Mode')


    def PlotCoupling(self, Scene, Col, Row, Combination):
        for (Mode0, Mode1) in Combination:
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.ITRList,
                          y        = self.GetCoupling( Mode0.ModeNumber, Mode1.ModeNumber ),
                          Fill     = False,
                          Legend   = f'{Mode0.ModeNumber} - {Mode1.ModeNumber}',
                          xLabel   = r'ITR',
                          yLabel   = r'Mode coupling')

            Scene.SetAxes(Col, Row, Equal=False, Legend=True, LegendTitle='Mode')


    def PlotAdiabatic(self, Scene, Col, Row, Combination):
        for (Mode0, Mode1) in Combination:
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.ITRList,
                          y        = self.GetAdiabatic( Mode0.ModeNumber, Mode1.ModeNumber ),
                          Fill     = False,
                          Legend   = f'{Mode0.ModeNumber} - {Mode1.ModeNumber}',
                          xLabel   = r'ITR',
                          yLabel   = r'Adiabatic criterion')

            Scene.SetAxes(Col, Row, Equal=False, Legend=True, yScale='log', yLimits=[1e-8, 1e-1], LegendTitle='Mode')


    def _PlotFields(self, iter=0):
        iter = ToList(iter)
        Scene = Scene2D(nCols=self.Size, nRows=len(iter), ColorBar=True, UnitSize=[5,5])

        for supermode in self.IterateSuperMode():
            for s, slice in enumerate(supermode.GetSlices(iter)):
                Field, x, y = slice.GetFullField()

                Scene.AddMesh(Row      = s,
                              Col      = supermode.ModeNumber,
                              x        = x,
                              y        = y,
                              Scalar   = Field,
                              ColorMap = FieldMap,
                              xLabel   = r'X-distance [$\mu$m]',
                              yLabel   = r'Y-distance [$\mu$m]' if supermode.ModeNumber==0 else "",
                              Title    = f'{supermode.Name} [ITR: {slice.ITR:.2f}]'
                              )

                Scene.SetAxes(supermode.ModeNumber, s, Equal=True)

        return Scene




    def PlotFields(self, iter=0):
        Scene =  self._PlotFields(iter=iter)

        Scene.Show()

    def GetCombination(self, Combination):
        if Combination is None:
            return tuple( combinations( self.SuperModes, 2 ) )
        else:
            Output = []
            for (c0, c1) in Combination:
                Output.append( ( self.SuperModes[c0], self.SuperModes[c1] ) )

    def _Plot(self, Input, iter=0, Combination=None):
        Input = ToList(Input)

        Combination = self.GetCombination(Combination)

        Scene = Scene2D(nCols=1, nRows=len(Input), ColorBar=False)

        i = 0
        if 'Index' in Input:
            self.PlotIndex(Scene, 0, i); i += 1

        if 'Beta' in Input:
            self.PlotBeta(Scene, 0, i); i += 1

        if 'Coupling' in Input:
            self.PlotCoupling(Scene, 0, i, Combination); i += 1

        if 'Adiabatic' in Input:
            self.PlotAdiabatic(Scene, 0, i, Combination); i += 1

        return Scene


    def Plot(self, Input, iter=0, Combination=None):
        Scene = self._Plot(Input, iter=iter, Combination=Combination)

        Scene.Show()


    def ExportPDF(self, Directory, iter=0, dpi=100, Combination=None):
        Scene0 = self._Plot(Input=['Index', 'Adiabatic'], iter=iter, Combination=Combination)
        Scene1 = self._PlotFields(iter=iter)

        dir = os.path.join(ZeroPath, Directory) + '.pdf'

        Multipage(dir, figs=[Scene0.Figure, Scene1.Figure], dpi=dpi)











# -
