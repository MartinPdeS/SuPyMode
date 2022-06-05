import logging
import os
import numpy               as np
import copy                as cp
from scipy.integrate       import solve_ivp
from scipy.interpolate     import interp1d
from mayavi                import mlab

from SuPyMode.Tools.BaseClass   import SetProperties, SetPlottings
from SuPyMode.Plotting.Plots      import Scene2D

Mlogger = logging.getLogger(__name__)


class ReprBase:
    Description = ''
    ReprVar     = []
    HLine       = "\n" + "--"*20 + "\n"
    Methods     = []

    def __repr__(self):
        String = []
        String.append(f"{self.Description}")
        String.append(self.HLine)
        String.append("Attributes:")

        for element in self.ReprVar:
            String += f"\n\t {element:20s}:\t {getattr(self, element)}"

        String.append("\n\nMethods:\n\t")
        String.append( f"\n\t".join(self.Methods) )
        String.append(self.HLine)
        return "".join(String)


class SuperSet(SetProperties, SetPlottings, ReprBase):

    Description = 'SuperSet class'

    ReprVar     = ["ParentSolver",
                   "Size",
                   "Geometry"]

    Methods     = [ "Propagate"]

    def __init__(self, ParentSolver):
        self.ParentSolver   = ParentSolver
        self.SuperModes     = []
        self._NextMode      = 0
        self._Matrix = None


    def ComputePropagationMatrix(self):
        M = np.zeros([self.Size, self.Size, len(self.ITRList)])
        for mode in self.SuperModes:
            M[mode.ModeNumber, mode.ModeNumber, :] = mode.Betas

        return M

    @property
    def Matrix(self):
        if self._Matrix is None:
            self._Matrix = self.ComputePropagationMatrix()
        return self._Matrix

    @property
    def NextMode(self):
        self._NextMode += 1
        return self._NextMode - 1

    def IterateSuperMode(self):
        for n, supermode in enumerate(self.SuperModes):
            yield supermode


    def ComputeM(self, CouplingFactor):
        shape = self.Beta.shape
        M     = np.zeros( [shape[0], shape[1], shape[1]] )
        for iter in range(shape[0]):
            beta = self.Beta[iter]
            M[iter] = CouplingFactor[iter] * self.Coupling[iter] + beta * np.identity(shape[1])

        return M


    def ComputeCouplingFactor(self, Length):
        dx =  Length/(self.Geometry.ITRList.size)

        dITR = np.gradient(np.log(self.Geometry.ITRList), 1)

        factor = dITR/dx

        return factor/10


    def Propagate(self, Amplitude=[1,1, 0, 0, 0], Length=1000):
        Amplitude = np.asarray(Amplitude)

        Distance = np.linspace(0, Length, self.ITRList.size)

        #Factor = self.ComputeCouplingFactor(Length)

        #M = self.ComputeM(CouplingFactor=Factor)



        Minterp = interp1d(Distance, self.Matrix, axis=-1)

        def foo(t, y):
            return 1j * Minterp(t).dot(y)

        sol = solve_ivp(foo,
                        y0       = Amplitude.astype(complex),
                        t_span   = [0, Length],
                        method   = 'RK45')

        return sol.y


    def Propagate_(self, Amplitude, Length, **kwargs):
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
                        **kwargs)

        return sol


    def PlotPropagation(self, Modes):
        Mode0 = self[0]
        Mode1 = self[1]

        Field = Mode0[0].GetFullField(Axes=False) + Mode1[0].GetFullField(Axes=False)

        surface = mlab.surf( np.abs( Field ) , warp_scale="auto" )

        @mlab.animate(delay=100)
        def anim_loc():
            for n, _ in Mode0.IterateSlice():
                Field = Mode0[n].GetFullField(Axes=False) + Mode1[n].GetFullField(Axes=False)

                surface.mlab_source.scalars = np.abs( Field )

                yield

        anim_loc()
        mlab.show


    @property
    def Size(self):
        return len(self.SuperModes)

    @property
    def Geometry(self):
        return self.ParentSolver.Geometry

    @property
    def ITRList(self):
        return self.ParentSolver.ITRList

    @property
    def Axes(self):
        return self.ParentSolver.Geometry.Axes


    def AppendSuperMode(self, CppSolver, BindingNumber):
        superMode = SuperMode(ParentSet=self, ModeNumber=self.NextMode, CppSolver=CppSolver, BindingNumber=BindingNumber )

        self.SuperModes.append( superMode )


    def __getitem__(self, N):
        return self.SuperModes[N]


    def __setitem__(self, N, val):
        self.SuperModes[N] = val











class SuperMode(ReprBase):
    Description = 'Supermode class'
    ReprVar     = ["ModeNumber",
                   "BindingNumber",
                   "ParentSet",
                   "LeftSymmetry",
                   "RightSymmetry",
                   "TopSymmetry",
                   "BottomSymmetry",
                   "Size"]

    Methods     = ["Fields",
                   "Index",
                   "Betas",
                   "PlotIndex",
                   "PlotBetas",
                   "PlotPropagation"]

    def __init__(self, ParentSet, ModeNumber, CppSolver,  BindingNumber):
        self.Binded         = CppSolver.GetMode(ModeNumber)
        self.ModeNumber     = ModeNumber

        self.CppSolver      = CppSolver
        self.ParentSet      = ParentSet

        self._Fields        = None
        self._FullFields    = None
        self._FullxAxis     = None
        self._FullxAxis     = None
        self._Index         = None
        self._Betas         = None


    @property
    def BindingNumber(self):
        return self.Binded.BindingNumber


    @property
    def Fields(self):
        if self._Fields is None:
            self._Fields = self.Binded.GetFields()
        return self._Fields

    @property
    def FullFields(self):
        if self._FullFields is None:
            self.ComputeFullFields()
        return self._FullFields


    @property
    def FullxAxis(self):
        if self._FullxAxis is None:
            self.ComputeFullFields()
        return self._FullxAxis


    @property
    def FullyAxis(self):
        if self._FullyAxis is None:
            self.ComputeFullFields()
        return self._FullyAxis


    @property
    def Index(self):
        if self._Index is None:
            self._Index = self.Binded.GetIndex()
        return self._Index


    @property
    def Betas(self):
        if self._Betas is None:
            self._Betas = self.Binded.GetBetas()
        return self._Betas


    def PlotIndex(self):
        Scene = Scene2D(nCols=1, nRows=1, ColorBar=False)

        Scene.AddLine(Row      = 0,
                      Col      = 0,
                      x        = self.ITRList,
                      y        = self.Index,
                      Fill     = False,
                      Legend   = self.ModeNumber,
                      xLabel   = r'ITR',
                      yLabel   = r'Effective refractive index n$_{eff}$')

        Scene.SetAxes(0, 0, Equal=False, Legend=True, yLimits=[self.Geometry.MinIndex/1.005, self.Geometry.MaxIndex], LegendTitle='Mode')
        Scene.Show()


    def PlotBetas(self):
        Scene = Scene2D(nCols=1, nRows=1, ColorBar=False)

        Scene.AddLine(Row      = 0,
                      Col      = 0,
                      x        = self.ITRList,
                      y        = self.Betas,
                      Fill     = False,
                      Legend   = self.ModeNumber,
                      xLabel   = r'ITR',
                      yLabel   = r'Propagation constante $\beta$')

        Scene.SetAxes(0, 0, Equal=False, Legend=True, LegendTitle='Mode')
        Scene.Show()


    def PlotFields(self, Slice: list):
        Scene = Scene2D(nCols=len(Slice), nRows=1, ColorBar=False, UnitSize=[3,3])

        for n, slice in enumerate(Slice):
            Scene.AddMesh(Row      = 0,
                          Col      = n,
                          x        = self.FullxAxis,
                          y        = self.FullyAxis,
                          Scalar   = self.FullFields[slice].T,
                          xLabel   = r'X-Direction [$\mu m$]',
                          yLabel   = r'Y-direction [$\mu m$]')

            Scene.SetAxes(Row=0, Col=n, Equal=True, Legend=False)
        Scene.Show()



    @property
    def LeftSymmetry(self):
        return self.Binded.LeftSymmetry

    @property
    def RightSymmetry(self):
        return self.Binded.RightSymmetry

    @property
    def TopSymmetry(self):
        return self.Binded.TopSymmetry

    @property
    def BottomSymmetry(self):
        return self.Binded.BottomSymmetry


    def IterateSlice(self):
        for n, slice in enumerate(self.Slice):
            yield n, slice

    @property
    def Size(self):
        return len(self.ParentSet.ITRList)

    @property
    def Geometry(self):
        return self.ParentSet.Geometry


    @property
    def ITRList(self):
        return self.ParentSet.ITRList


    @property
    def Axes(self):
        return self.ParentSet.Axes


    def ApplySymmetry(self, Field, xAxis, yAxis, axis):
        Slice = [slice(None), slice(None), slice(None, None, -1)]
        Output = np.concatenate((Field[Slice], Field), axis=axis)


    def ExtendAxis(self, Axis: np.ndarray, sign: str):
        d     = Axis[1] - Axis[0]

        if sign == "Plus":
            start = Axis[-1] + d
            Next  = np.arange(0, Axis.size) * d + start
            Next  = [Axis, Next]

        if sign == "Minus":
            stop = Axis[0]
            Next  = np.arange(-Axis.size, 0) * d + stop
            Next  = [Next, Axis]

        return np.concatenate(Next)



    def ComputeFullFields(self):
        self._FullxAxis = self.Axes.X
        self._FullyAxis = self.Axes.Y
        self._FullFields = self.Fields

        if self.BottomSymmetry == 1:
            self._FullFields = np.concatenate((self._FullFields[:, :, ::-1], self._FullFields), axis=2)
            self._FullyAxis = self.ExtendAxis(Axis=self._FullyAxis, sign="Minus")

        elif self.BottomSymmetry == -1:
            self._FullFields = np.concatenate((-self._FullFields[:, :, ::-1], self._FullFields), axis=2)
            self._FullyAxis = self.ExtendAxis(Axis=self._FullyAxis, sign="Minus")

        if self.TopSymmetry == 1:
            self._FullFields = np.concatenate((self._FullFields, self._FullFields[:, :, ::-1]), axis=2)
            self._FullyAxis = self.ExtendAxis(Axis=self._FullyAxis, sign="Plus")

        elif self.TopSymmetry == -1:
            self._FullFields = np.concatenate((self._FullFields, -self._FullFields[:, :, ::-1]), axis=2)
            self._FullyAxis = self.ExtendAxis(Axis=self._FullyAxis, sign="Plus")


        if self.RightSymmetry == 1:
            self._FullFields = np.concatenate((self._FullFields[...], self._FullFields[:, ::-1, :]), axis=1) #Here
            self._FullxAxis = self.ExtendAxis(Axis=self._FullxAxis, sign="Plus")

        elif self.RightSymmetry == -1:
            self._FullFields = np.concatenate((self._FullFields[...], -self._FullFields[:, ::-1, :]), axis=1) #Here
            self._FullxAxis = self.ExtendAxis(Axis=self._FullxAxis, sign="Plus")

        if self.LeftSymmetry == 1:
            self._FullFields = np.concatenate((self._FullFields[:, ::-1, :], self._FullFields[...]), axis=1) #Here
            self._FullxAxis = self.ExtendAxis(Axis=self._FullxAxis, sign="Minus")

        elif self.LeftSymmetry == -1:
            self._FullFields = np.concatenate((-self._FullFields[:, ::-1, :], self._FullFields[...]), axis=1) #Here
            self._FullxAxis = self.ExtendAxis(Axis=self._FullxAxis, sign="Minus")



    def CompareSymmetries(self, Other):
        assert isinstance(Other, SuperMode), "Can only compare SuperMode instance with another"
        return np.all(self.Symmetries == Other.Symmetries)


    def AppendSlice(self, SliceNumber):
        self.Slice.append( SetSlice(ParentMode=self, SliceNumber=SliceNumber) )

    def GetSlices(self, SlicesNumber: list):
        return [self[slice] for slice in SlicesNumber]


    def __getitem__(self, N):
        return self.Slice[N]


    def __setitem__(self, N, val):
        self.Slice[N] = val


    def PlotPropagation(self):
        surface = mlab.surf( np.abs( self.Fields[0]), warp_scale="auto")

        @mlab.animate(delay=100)
        def anim_loc():
            for field in self.Fields:
                surface.mlab_source.scalars = np.abs( field )

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





# -
