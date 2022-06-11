import logging
import os
import numpy               as np
import copy                as cp
from scipy.integrate       import solve_ivp
from scipy.interpolate     import interp1d
from mayavi                import mlab

from SuPyMode.Tools.BaseClass   import SetProperties, SetPlottings
from SuPyMode.Plotting.Plots      import Scene, Axis, Line, Mesh
from SuPyMode.Plotting.PlotsUtils  import FieldMap, MidPointNorm
from SuPyMode.Tools.Directories import RootPath

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



class SuperPosition(ReprBase):
    Description = 'Mode superposition class'
    ReprVar     = ["Amplitudes"]

    def __init__(self, SuperSet, InitialAmplitudes: list):
        self.SuperSet   = SuperSet
        self.InitialAmplitudes = np.asarray(InitialAmplitudes).astype(complex)
        self._CouplerLength    = None
        self._Amplitudes       = None


    def ComputeAmpltiudes(self, rTol=1e-8, aTol=1e-7, MaxStep=np.inf):
        self.MatrixInterp = interp1d(self.Distance, self.SuperSet.Matrix, axis=-1)

        def foo(t, y):
            return 1j * self.MatrixInterp(t).dot(y)

        sol = solve_ivp(foo,
                        y0       = self.InitialAmplitudes,
                        t_span   = [0, self._CouplerLength],
                        method   = 'RK45',
                        rtol     = rTol,
                        atol     = aTol,
                        max_step = MaxStep)

        self.Amplitudes = sol.y
        self.Distances  = sol.t
        self.AmplitudeInterpolation = interp1d(self.Distances, self.Amplitudes, axis=-1)

    @property
    def ITRList(self):
        return self.SuperSet.ITRList

    @property
    def CouplerLength(self):
        assert self._CouplerLength is not None, "CouplerLength attribute has to be defined before computing propagation."
        return self._CouplerLength

    @CouplerLength.setter
    def CouplerLength(self, Value: float):
        self._CouplerLength = Value
        self.Distance = np.linspace(0, self._CouplerLength, self.ITRList.size)


    @property
    def Amplitudes(self):
        if self._Amplitudes is None:
            self.ComputeAmpltiudes()

        return self._Amplitudes


    @Amplitudes.setter
    def Amplitudes(self, Value):
        self._Amplitudes = Value


    def PlotAmplitudes(self):
        Scene = Scene2D(nCols=1, nRows=1, ColorBar=False)
        A = self.InitialAmplitudes.dot(self.Amplitudes)
        z = self.Distances
        Scene.AddLine(Row      = 0,
                      Col      = 0,
                      x        = z,
                      y        = A.real,
                      Fill     = False,
                      Legend   = r"$\Re${A}",
                      xLabel   = r'Z-Distance [$\mu m$]',
                      yLabel   = r'Mode complex ampltiude [normalized]')

        Scene.AddLine(Row      = 0,
                      Col      = 0,
                      x        = z,
                      y        = np.abs(A),
                      Fill     = False,
                      Legend   = r"|A|",
                      xLabel   = r'Z-Distance [$\mu m$]',
                      yLabel   = r'Mode complex ampltiude [normalized]')


        Scene.SetAxes(0, 0, Equal=False, Legend=True)
        Scene.Show()

    def PlotField(self, Slice):
        if self._Amplitudes is None: self.ComputeAmpltiudes()

        Scene = Scene2D(nCols=1, nRows=1, ColorBar=False)

        Field = self.SuperSet[0].FullFields[0]*0.

        for a, mode in zip(self.InitialAmplitudes, self.SuperSet.SuperModes):
            field = np.real(a)*mode.FullFields[Slice]
            Field += field


        Scene.AddMesh(Row      = 0,
                      Col      = 0,
                      x        = self.SuperSet[0].FullxAxis,
                      y        = self.SuperSet[0].FullyAxis,
                      Scalar   = Field.T.real,
                      xLabel   = r'X-Direction [$\mu m$]',
                      yLabel   = r'Y-direction [$\mu m$]')

        Scene.SetAxes(Row=0, Col=0, Equal=True, Legend=False)
        Scene.Show()


    def PlotPropagation(self):
        if self._Amplitudes is None: self.ComputeAmpltiudes()

        y = self.AmplitudeInterpolation(self.Distance)

        z = self.Distance

        Field = self.SuperSet[0].FullFields.astype(complex)*0.

        for mode, _ in enumerate(self.InitialAmplitudes):
            a = y[mode].astype(complex)
            field = self.SuperSet[mode].FullFields.astype(complex)
            Field += np.einsum('i, ijk->ijk', a, field)

        surface = mlab.surf( np.abs( Field[0] ) , warp_scale="auto" )

        @mlab.animate(delay=100)
        def anim_loc():
            for n, _ in enumerate(self.Distance):
                surface.mlab_source.scalars = np.abs(np.abs( Field[n] ) )

                yield

        anim_loc()
        mlab.show()





class SuperSet(SetProperties, SetPlottings, ReprBase):

    Description = 'SuperSet class'

    ReprVar     = ["ParentSolver",
                   "Size",
                   "Geometry"]

    Methods     = ["GetSuperposition", "Matrix"]

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

        return dITR/dx


    def GetSuperposition(self, Amplitudes):
        return SuperPosition(SuperSet=self, InitialAmplitudes=Amplitudes)


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


    def AppendSuperMode(self, CppSolver, BindingNumber, SolverNumber):
        superMode = SuperMode(ParentSet=self, ModeNumber=self.NextMode, CppSolver=CppSolver, BindingNumber=BindingNumber, SolverNumber=SolverNumber )

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

    def __init__(self, ParentSet, ModeNumber, CppSolver,  BindingNumber, SolverNumber):
        self.Binded         = CppSolver.GetMode(BindingNumber)
        self.ModeNumber     = ModeNumber
        self.SolverNumber   = SolverNumber
        self.ID             = [SolverNumber, BindingNumber]
        self.Name           = f"Mode {self.ID[0]}:{self.ID[1]}"

        self.CppSolver      = CppSolver
        self.ParentSet      = ParentSet

        self._Fields        = None
        self._FullFields    = None
        self._FullxAxis     = None
        self._FullxAxis     = None
        self._Index         = None
        self._Betas         = None
        self._Adiabatic     = None
        self._Coupling      = None


    @property
    def BindingNumber(self):
        return self.Binded.BindingNumber


    @property
    def Adiabatic(self):
        if self._Adiabatic is None:
            self._Adiabatic = self.Binded.GetAdiabatic()

        return self._Adiabatic


    @property
    def Coupling(self):
        if self._Coupling is None:
            self._Coupling = self.Binded.GetCoupling()
        return self._Coupling


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


    def _PlotIndex(self, Ax):
        artist = Line(X=self.ITRList, Y=self.Index, Label=None, Fill=False)

        Ax.AddArtist(artist)


    def _PlotBetas(self, Ax):
        artist = Line(X=self.ITRList, Y=self.Betas, Label=self.Name, Fill=False)

        Ax.AddArtist(artist)


    def PlotIndex(self):
        Fig = Scene('SuPyMode Figure', UnitSize=(10,4))

        ax = Axis(Row    = 0,
                  Col    = 0,
                  xLabel = 'ITR',
                  yLabel = r'Effective refraction index',
                  Title  = None,
                  Grid   = True,
                  xScale = 'linear',
                  yScale = 'linear')

        self._PlotIndex(ax)

        Fig.AddAxes(ax)

        Fig.Show()


    def PlotBetas(self):
        Fig = Scene('SuPyMode Figure', UnitSize=(10,4))

        ax = Axis(Row    = 0,
                  Col    = 0,
                  xLabel = 'ITR',
                  yLabel = r'Propagation constante $\beta$',
                  Title  = None,
                  Grid   = True,
                  xScale = 'linear',
                  yScale = 'linear')

        self._PlotBetas(ax)

        Fig.AddAxes(ax)

        Fig.Show()


    def _PlotFields(self, Ax, slice):
        artist = Mesh(X           = self.FullxAxis,
                      Y           = self.FullyAxis,
                      Scalar      = self.FullFields[slice].T,
                      ColorMap    = FieldMap,
                      DiscretNorm = False,
                      )

        Ax.AddArtist(artist)


    def PlotFields(self, Slice: list):
        Fig = Scene('SuPyMode Figure', UnitSize=(10,4))

        for n, slice in enumerate(Slice):
            ax = Axis(Row    = 0,
                      Col    = n,
                      xLabel = r'X-Direction [$\mu m$]',
                      yLabel = r'Y-direction [$\mu m$]',
                      Title  = f'{self.Name}  [ITR: {self.ITRList[slice]:.2f}]',
                      Legend = False,
                      ColorBar=True,
                      ColorbarPosition = 'right',
                      Grid   = True,
                      Equal  = True,
                      xScale = 'linear',
                      yScale = 'linear')

            self._PlotFields(ax, slice)

            Fig.AddAxes(ax)

        Fig.Show()


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


    def GetSlice(self, Slice: int, Full: bool=True):
        if Full:
            return self.Betas[Slice], self.FullFields[Slice], self._FullxAxis, self._FullyAxis
        else:
            return self.Betas[Slice], self.Fields[Slice], self.Axis.X, self.Axis.Y

    def __getitem__(self, N):
        return self.Slice[N]


    def __setitem__(self, N, val):
        self.Slice[N] = val


    def GetArrangedFields(self):
        sign = np.sign( np.sum(self.FullFields[0]))
        FullFields = [sign*self.FullFields[0]]

        for field in self._FullFields:
            overlap = np.sum(field*FullFields[-1])
            if overlap > 0:
                FullFields.append(field/np.max(np.abs(field)))

            if overlap <= 0:
                FullFields.append(-field/np.max(np.abs(field)))

        return FullFields

    def PlotPropagation(self, SaveName=None):

        FullFields = self.GetArrangedFields()

        FileName = []

        factor = 5
        offset = 11

        fig = mlab.figure(size=(1000,700), bgcolor=(1,1,1), fgcolor=(0,0,0))

        surface = mlab.surf(FullFields[0]*factor + offset, colormap='coolwarm', warp_scale='4', representation='wireframe', line_width=6, opacity=0.9, transparent=True)

        mesh = self.Geometry.GetFullMesh(self.LeftSymmetry, self.RightSymmetry, self.TopSymmetry, self.BottomSymmetry)
        baseline = mlab.surf(mesh*0, color=(0,0,0), representation='wireframe', opacity=0.53)

        #mlab.contour_surf(mesh, color=(0,0,0), contours=[mesh.min(), 1.4, mesh.max()], line_width=6)

        mlab.axes( xlabel='x', ylabel='y', zlabel='z', color=(0,0,0), nb_labels=10, ranges=(0,40,0,40,0,20), y_axis_visibility=False )


        mlab.gcf().scene.parallel_projection = False
        mlab.view(elevation=70, distance=300)
        mlab.move(up=-6)

        #mlab.outline(baseline)


        import imageio

        @mlab.animate(delay=10)
        def anim_loc():
            for n, field in enumerate(FullFields):
                surface.mlab_source.scalars = field*factor + offset
                baseline.mlab_source.scalars = field*3


                if SaveName is not None:
                    FileName.append( f'{RootPath}/Animation/Animation_{n:03d}.png' )
                    mlab.savefig(filename=FileName[-1])

                yield

        anim_loc()
        mlab.show()


        if SaveName is not None:

            with imageio.get_writer(f'{RootPath}/Animation/{SaveName}.gif', mode='I', fps=50) as writer:
                for filename in FileName:
                    image = imageio.imread(filename)
                    writer.append_data(image)






# -
