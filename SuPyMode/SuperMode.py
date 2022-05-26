import logging
import numpy               as np
import copy                as cp
from scipy.integrate       import solve_ivp
from scipy.interpolate     import interp1d
from itertools             import combinations, combinations_with_replacement as Combinations
from mayavi                import mlab

from SuPyMode.Tools.Config      import *
from SuPyMode.Tools.utils       import RecomposeSymmetries, Enumerate
from SuPyMode.Tools.BaseClass   import SetProperties
from SuPyMode.Tools.utils       import Multipage, prePlot, ToList
from SuPyMode.Plotting.Plots       import Scene2D
from SuPyMode.Plotting.PlotsUtils  import FieldMap


Mlogger = logging.getLogger(__name__)

class SuperMode(object):

    def __init__(self, Name, Geometry, ParentSet):
        self.Name       = Name
        self.Slice      = []
        self._Adiabatic = None
        self._Coupling  = None
        self.Geometry   = Geometry
        self.Parent     = ParentSet



    def Append(self, **kwargs):
        self.Slice.append(kwargs['Field'])


    def __getitem__(self, N):
        return self.Slice[N]


    def __setitem__(self, N, val):
        self.Slice[N] = val


    def PlotPropagation(self):
        image, _, _ = self.FullField(0)

        image = np.abs(image)

        mlab.surf(image, warp_scale="auto")

        @mlab.animate(delay=100)
        def anim_loc():
            for i in range(len(self.Field)):
                image, _, _ = self.FullField(i)
                mlab.surf( np.abs(image), warp_scale="auto")
                yield

        anim_loc()
        mlab.show()
        """
        import os
        fps = 20
        prefix = 'ani'
        ext = '.png'

        import subprocess
        animate_plots(base_directory='yolo', fname_prefix='dasda')"""


    def FullField(self, iter):
        Field      = self.Slice[iter]

        Field, xAxes, yAxes = RecomposeSymmetries(Input      = Field,
                                                  Symmetries = ParentSet.Symmetries,
                                                  Axes       = ParentSet.Geometry.Axes)

        return Field, xAxes, yAxes


    def __copy__(self):
        to_be_copied = ['Slice']

        copy_ = SuperSet( self.Name, self.Geometry)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_


    def __deepcopy__(self, memo):
        to_be_copied = ['Slice']

        copy_ = SuperMode( self.Name, self.Geometry)

        for attr in self.__dict__:

            if attr in to_be_copied:
                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])

            else:
                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_



class _SuperMode(object):

    def __init__(self, Amplitude, Name, ParentSet):
        self.InitAmplitude = np.asarray(Amplitude).astype(complex)
        self.Name          = Name
        self.Parent        = ParentSet
        self.Amplitudes    = None
        self.Amplitude     = None


    def ComputeCouplingFactor(self, Profile, Distance):
        dx =  np.gradient(Distance, 1)

        dITR = np.gradient(np.log(Profile), 1)

        factor = dITR/dx

        return factor


    def Propagate(self, Profile, Distance, **kwargs):

        Profile = np.asarray(Profile)

        Factor = self.ComputeCouplingFactor(Profile, Distance)

        M = self.Parent.ComputeM(CouplingFactor=Factor)

        Minterp = interp1d(Distance, M, axis=0)

        def foo(t, y):
            return 1j * Minterp(t).dot(y)

        Sol = solve_ivp(foo,
                        y0       = self.InitAmplitude,
                        t_span   = [Distance[0], Distance[-1]],
                        method   = 'RK45',
                        vectorized=True,
                        **kwargs
                        )

        self.Z          = Sol.t
        self.Amplitudes = Sol.y



    def PlotPropagation(self):
        image, _, _ = self.FullField(0)

        image = np.abs(image)

        mlab.surf(image, warp_scale="auto")

        @mlab.animate(delay=100)
        def anim_loc():
            for i in range(len(self.Field)):
                image, _, _ = self.FullField(i)
                mlab.surf( np.abs(image), warp_scale="auto")
                yield

        anim_loc()
        mlab.show()
        """
        import os
        fps = 20
        prefix = 'ani'
        ext = '.png'

        import subprocess
        animate_plots(base_directory='yolo', fname_prefix='dasda')"""


    def FullField(self, iter):
        Field      = self.Slice[iter]

        Field, xAxes, yAxes = RecomposeSymmetries(Input      = Field,
                                                  Symmetries = ParentSet.Symmetries,
                                                  Axes       = ParentSet.Geometry.Axes)

        return Field, xAxes, yAxes


    def __copy__(self):
        to_be_copied = ['Slice']

        copy_ = SuperSet( self.Name, self.Geometry)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_


    def __deepcopy__(self, memo):
        to_be_copied = ['Slice']

        copy_ = SuperMode( self.Name, self.Geometry)

        for attr in self.__dict__:

            if attr in to_be_copied:
                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])

            else:
                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_



class SuperSet(SetProperties):
    def __init__(self, NSolutions, Geometry, debug='INFO'):
        Mlogger.setLevel(getattr(logging, debug))

        self.sMode    = NSolutions
        self.SuperModes    = []
        self._Coupling     = None
        self._Adiabatic    = None
        self._Index        = None
        self._Beta         = None
        self.Geometry      = Geometry
        self._M            = None

        self.Init()

        self.combinations = tuple(combinations( np.arange(NSolutions), 2 ) )
        self.Combinations = tuple(Combinations( np.arange(NSolutions), 2 ) )


    def PlotIndex(self, Scene, Col, Row):
        for i in range(self.sMode):
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.Geometry.ITRList,
                          y        = self.Index[:,i],
                          Fill     = False,
                          Legend   = f'Mode: {i}',
                          xLabel   = r'ITR',
                          yLabel   = r'Effective index n_{eff}',
                          )

            Scene.SetAxes(0, 0, Equal=False, Legend=True)


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

            Scene.SetAxes(0, 0, Equal=False, Legend=True)


    def PlotCoupling(self, Scene, Col, Row):
        for n, (i,j) in enumerate( self.Combination ):
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.Geometry.ITRList,
                          y        = self.Coupling[:,i,j],
                          Fill     = False,
                          Legend   = f'Mode: {i}-{j}',
                          xLabel   = r'ITR',
                          yLabel   = r'Mode coupling')

            Scene.SetAxes(0, 0, Equal=False, Legend=True)


    def PlotAdiabatic(self, Scene, Col, Row):
        for n, (i,j) in enumerate( self.Combination ):
            Scene.AddLine(Row      = Row,
                          Col      = Col,
                          x        = self.Geometry.ITRList,
                          y        = self.Adiabatic[:,i,j],
                          Fill     = False,
                          Legend   = f'Mode: {i}-{j}',
                          xLabel   = r'ITR',
                          yLabel   = r'Adiabatic criterion')

            Scene.SetAxes(0, 0, Equal=False, Legend=True)


    def PlotFields(self, iter=0):
        Scene = Scene2D(nCols=self.sMode, nRows=1, ColorBar=False, UnitSize=[5,5])

        for mode in range(self.sMode):
            Field, x, y = self[mode][iter].GetField()

            Scene.AddMesh(Row      = 0,
                          Col      = mode,
                          x        = x,
                          y        = y,
                          Scalar   = Field,
                          ColorMap = FieldMap,
                          xLabel   = r'X-distance [$\mu$m]',
                          yLabel   = r'Y-distance [$\mu$m]',
                          Title    = f'Mode: {mode} [ITR: {self.Geometry.ITRList[iter]}]'
                          )

            Scene.SetAxes(mode, 0, Equal=True)

        Scene.Show()



    def Plot(self, Input, iter=0, PlotKwarg=None, Combination=None):
        Input = ToList(Input)

        if Combination is None:
            self.Combination = tuple(combinations( np.arange(self.sMode), 2 ) )
        else:
            self.Combination = Combination

        Scene = Scene2D(nCols=1, nRows=len(Input), ColorBar=False)

        i = 0
        if 'Index' in Input:
            self.PlotIndex(Scene, 0, i); i += 1

        if 'Beta' in Input:
            self.PlotBeta(Scene, 0, i); i += 1

        if 'Coupling' in Input:
            self.PlotCoupling(Scene, 0, i); i += 1

        if 'Adiabatic' in Input:
            self.PlotAdiabatic(Scene, 0, i); i += 1

        Scene.Show()



    def SaveFig(self, Input, Directory, iter=0, dpi=100, PlotKwarg=None, Combination=None):

        if Combination is None:
            self.Combination = tuple(combinations( np.arange(self.sMode), 2 ) )
        else:
            self.Combination = Combination

        figures = self.GenFigures(Input, iter, PlotKwarg)

        dir = os.path.join(ZeroPath, Directory) + '.pdf'

        Multipage(dir, figs=figures, dpi=dpi)


    def ComputeCoupling(self):
        return self.CppSolver.ComputingCoupling()


    def Init(self):
        for solution in range(self.sMode):
            supermode = SuperMode(Name = f"Mode {solution}", Geometry = self.Geometry, ParentSet=self)
            supermode.number = solution
            self.SuperModes.append(supermode)


    def ComputeM(self, CouplingFactor):
        shape = self.Beta.shape
        M     = np.zeros( [shape[0], shape[1], shape[1]] )
        for iter in range(shape[0]):
            beta = self.Beta[iter]
            M[iter] = CouplingFactor[iter] * self.Coupling[iter] + beta * np.identity(shape[1])

        return M


    def GetMode(self, Amplitude, Name = 1):
        return _SuperMode(Amplitude=Amplitude, Name=Name, ParentSet=self)


    def ComputeCouplingFactor(self, Length):
        dx =  Length/(self.Geometry.ITRList.size)

        dITR = np.gradient(np.log(self.Geometry.ITRList), 1)

        factor = dITR/dx

        return factor/10


    def Propagate(self, Amplitude, Length, **kwargs):
        Amplitude = np.asarray(Amplitude)

        Distance = np.linspace(0, Length, self.Geometry.ITRList.size)

        Factor = self.ComputeCouplingFactor(Length)


        M = self.ComputeM(CouplingFactor=Factor)

        Minterp = interp1d(Distance, M, axis=0)

        def foo(t, y):
            return 1j * Minterp(t).dot(y)

        sol = solve_ivp(foo,
                        y0       = Amplitude.astype(complex),
                        t_span   = [0, Length],
                        method   = 'RK45',
                        **kwargs
                        )

        return sol


    def __getitem__(self, N):
        return self.SuperModes[N]


    def __setitem__(self, N, val):
        self.SuperModes[N] = val


    def __copy__(self):
        to_be_copied = ['SuperModes']

        copy_ = SuperSet(self.sMode, self.Geometry)

        for attr in self.__dict__:

            if attr in to_be_copied:
                copy_.__dict__[attr] =  cp.deepcopy( self.__dict__[attr] )
            else:
                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_





class SetSlice(np.ndarray):

    def __new__(cls, Field, Index, Beta, ParentSet):
        self      = Field.view(SetSlice)

        return self


    def __init__(self, Field, Index, Beta, ParentSet):

        self.Field  = Field
        self.Index  = Index
        self.Beta   = Beta
        self.ParentSet = ParentSet


    def __array_finalize__(self, viewed):
        pass


    def __pow__(self, other):
        assert isinstance(other, SetSlice), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    def Overlap(self, other):
        assert isinstance(other, SetSlice), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    @property
    def Norm(self):
        return (self.Field**2).sum()

    def Normalize(self):
        norm = self.Norm
        self.Field *= 1 / norm
        self.Index *= 1 / norm

    def GetField(self, Symmetries=True):
        if Symmetries:
            return RecomposeSymmetries(self.Field, self.ParentSet.Symmetries, self.ParentSet.Geometry.Axes)
        else:
            self.Field, self.ParentSet.Geometry.Axes


    def __plot__(self, ax, title=None):

        Field, xaxis, yaxis = RecomposeSymmetries(self.Field, self.ParentSet.Symmetries, self.ParentSet.Geometry.Axes)


        ax.pcolormesh(yaxis, xaxis, Field.T, shading='auto', cmap='bwr')

        ax.set_ylabel(r'Y-distance [$\mu$m]', fontsize=6)

        ax.set_xlabel(r'X-distance [$\mu$m]', fontsize=6)

        ax.set_aspect('equal')
        if title:
            ax.set_title(title, fontsize=8)


    def __copy__(self):
        to_be_copied = ['Field', 'Index', 'Axes', 'Beta']

        copy_ = SetSlice(self, self.Axes, self.Index, self.Beta)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_



    def __deepcopy__(self, memo):
        to_be_copied = ['Field', 'Index', 'Axes', 'Beta']

        copy_ = SetSlice(self, self.Axes, self.Index, self.Beta)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_


















# -
