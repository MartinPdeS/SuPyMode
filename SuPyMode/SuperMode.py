import logging
import numpy               as np
from scipy.integrate       import solve_ivp
from scipy.interpolate     import interp1d
from mayavi                import mlab


from SuPyMode.Plotting.Plots      import Scene, Axis, Line, Mesh, ColorBar
from SuPyMode.Plotting.PlotsUtils import FieldMap
from SuPyMode.Tools.Directories   import RootPath
from SuPyMode.Tools.BaseClass     import ReprBase, ExtendField

Mlogger = logging.getLogger(__name__)


class SuperMode(ReprBase, ExtendField):
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
        self.Name           = f"Mode {SolverNumber}:{BindingNumber}"

        self.CppSolver      = CppSolver
        self.ParentSet      = ParentSet

        self._Fields        = None
        self._Index         = None
        self._Betas         = None
        self._Adiabatic     = None
        self._Coupling      = None

    @property
    def FullFields(self):
        if self._FullFields is None:
            self.ComputeFullFields()
        return self._FullFields


    @property
    def FullxAxis(self):
        if self._FullxAxis is None:
            self._FullxAxis, self._FullyAxis = self.GetFullAxis(self.Axes.X, self.Axes.Y)
        return self._FullxAxis

    @property
    def FullyAxis(self):
        if self._FullyAxis is None:
            self._FullxAxis, self._FullyAxis = self.GetFullAxis(self.Axes.X, self.Axes.Y)
        return self._FullyAxis

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
                      )

        Ax.AddArtist(artist)


    def PlotFields(self, Slice: list):
        Fig = Scene('SuPyMode Figure', UnitSize=(10,4))

        Colorbar = ColorBar(Discreet=False, Position='right')

        for n, slice in enumerate(Slice):
            ax = Axis(Row      = 0,
                      Col      = n,
                      xLabel   = r'X-Direction [$\mu m$]',
                      yLabel   = r'Y-direction [$\mu m$]',
                      Title    = f'{self.Name}  [ITR: {self.ITRList[slice]:.2f}]',
                      Legend   = False,
                      Colorbar = Colorbar,
                      Grid     = True,
                      Equal    = True,
                      xScale   = 'linear',
                      yScale   = 'linear')

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
