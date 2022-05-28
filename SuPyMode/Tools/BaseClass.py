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
    def ITR(self):
        return self.Geometry.ITRList

    @property
    def Symmetries(self):
        return self.Geometry.Symmetries

    @property
    def CppSolver(self):
        return self.Parent.CppSolver


    @property
    def LeftSymmetry(self):
        return self.Parent.LeftSymmetry


    @property
    def RightSymmetry(self):
        return self.Parent.RightSymmetry


    @property
    def TopSymmetry(self):
        return self.Parent.TopSymmetry


    @property
    def BottomSymmetry(self):
        return self.Parent.BottomSymmetry


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
        for i in range(self.sMode):
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.Geometry.ITRList,
                          y        = self.Index[:,i],
                          Fill     = False,
                          Legend   = f'Mode: {i}',
                          xLabel   = r'ITR',
                          yLabel   = r'Effective index n$_{eff}$',
                          )

            Scene.SetAxes(Col, Row, Equal=False, Legend=True, yLimits=[self.Parent.Geometry.MinIndex/1.005, self.Parent.Geometry.MaxIndex])


    def PlotBeta(self, Scene, Col, Row):
        for i in range(self.sMode):
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.Geometry.ITRList,
                          y        = self.Beta[:,i],
                          Fill     = False,
                          Legend   = f'Mode: {i}',
                          xLabel   = r'ITR',
                          yLabel   = r'Propagation constante $\beta$')

            Scene.SetAxes(Col, Row, Equal=False, Legend=True)


    def PlotCoupling(self, Scene, Col, Row, Combination):
        for n, (i,j) in enumerate( Combination ):
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.Geometry.ITRList,
                          y        = self.Coupling[:,i,j],
                          Fill     = False,
                          Legend   = f'Mode: {i}-{j}',
                          xLabel   = r'ITR',
                          yLabel   = r'Mode coupling')

            Scene.SetAxes(Col, Row, Equal=False, Legend=True)


    def PlotAdiabatic(self, Scene, Col, Row, Combination):
        for n, (i,j) in enumerate( Combination ):
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.Geometry.ITRList,
                          y        = self.Adiabatic[:,i,j],
                          Fill     = False,
                          Legend   = f'Mode: {i}-{j}',
                          xLabel   = r'ITR',
                          yLabel   = r'Adiabatic criterion')

            Scene.SetAxes(Col, Row, Equal=False, Legend=True, yScale='log', yLimits=[1e-8, 1e-1])


    def _PlotFields(self, iter=0):
        iter = ToList(iter)
        Scene = Scene2D(nCols=self.sMode, nRows=len(iter), ColorBar=True, UnitSize=[5,5])


        for m, supermode in self.IterateSuperMode():
            for s, slice in enumerate(iter):
                Field, x, y = supermode[s].GetField()
                Scene.AddMesh(Row      = s,
                              Col      = m,
                              x        = x,
                              y        = y,
                              Scalar   = Field,
                              ColorMap = FieldMap,
                              xLabel   = r'X-distance [$\mu$m]',
                              yLabel   = r'Y-distance [$\mu$m]' if m==0 else "",
                              Title    = f'{supermode.Name} [ITR: {self.Geometry.ITRList[slice]:.2f}]'
                              )

                Scene.SetAxes(m, s, Equal=True)




        return Scene




    def PlotFields(self, iter=0):
        Scene =  self._PlotFields(iter=iter)

        Scene.Show()


    def _Plot(self, Input, iter=0, Combination=None):
        Input = ToList(Input)

        if Combination is None: Combination = tuple(combinations( np.arange(self.sMode), 2 ) )

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
