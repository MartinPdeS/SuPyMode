
""" standard imports """
import os
import numpy as np
import logging

from matplotlib.path             import Path
from shapely.geometry            import Point
from shapely                     import affinity
from scipy.ndimage.filters       import gaussian_filter
import pickle

""" package imports """
from SuPyMode.Tools.Directories       import RootPath
from SuPyMode.Tools.Special           import gradientO4
from SuPyMode.Tools.utils             import ToList, Axes
from SuPyMode.Plotting.Plots          import Scene2D

Mlogger = logging.getLogger(__name__)

class Gradient:
    def __init__(self, Center, Nin, Nout, Rout):
        if Nin == Nout:
            Nin = Nout + 1e-9
        self.Nin    = Nin
        self.Nout   = Nout
        self.Rout   = Rout
        if isinstance(Center, Point):
            Center = Center.xy
        self.Center = Center

    def Evaluate(self, Xmesh, Ymesh):
        A = -(self.Rout * self.Nout) / (self.Nin-self.Nout)
        B = - A * self.Nin
        rho = np.hypot( Xmesh-self.Center[0], Ymesh-self.Center[1] )
        return B/(rho-A)


class Geometry(object):
    """ Class represent the refractive index (RI) geometrique profile which
    can be used to retrieve the supermodes.

    Parameters
    ----------
    Objects : type
        Geometrique object representing the element of the RI profile.
    Xbound : :class:`list`
        X-dimension boundary for the rasterized profile.
    Ybound : :class:`list`
        Y-dimension boundary for the rasterized profile.
    Nx : :class:`int`
        Number of point for X dimensions discretization.
    Ny : :class:`int`
        Number of point for Y dimensions discretization.
    """

    def __init__(self, Clad, Objects, Xbound, Ybound, Nx, Ny, debug='INFO'):

        Mlogger.setLevel(getattr(logging, debug))

        self.Clad       = Clad

        self.Objects    = ToList(Objects)

        #self.ScaleObject()

        self.Boundaries = [Xbound, Ybound]

        self.Shape      = [Nx, Ny]

        self.Axes  = Axes( {'wavelength': 1.0,
                            'Xbound'    : Xbound,
                            'Ybound'    : Ybound,
                            'Nx'        : Nx,
                            'Ny'        : Ny } )


    def ScaleObject(self):
        for object in self.Objects:
            Factor = object.Radius/self.Clad.Radius
            object.Scale(Factor=Factor)



    def Rotate(self, angle):
        for object in self.Objects:
            object.Object = affinity.rotate(object.Object, 28, (0,0))


    def rasterize_polygone(self, polygone):
        """ The method rasterize individual polygone object.

        Parameters
        ----------
        polygone : :class:`geo`
            The polygone to be rasterize.

        Returns
        -------
        :class:`np.ndarray`
            The rasterized object.

        """

        obj_ext = Path(list( polygone.Object.exterior.coords))

        obj_ext = obj_ext.contains_points(self.coords).reshape(self.Shape)

        if polygone.hole is not None:
            obj_int = Path(list( polygone.hole.exterior.coords))

            obj_int = np.logical_not( obj_int.contains_points(self.coords).reshape(self.Shape) )

            obj_ext = np.logical_and( obj_ext, obj_int )

        obj_ext = obj_ext.astype(float)

        polygone.raster = obj_ext



    def add_object_to_mesh(self, polygone):
        """ Method add a rasterized object to the pre-existing profile mesh
        with its provided refractive index.

        Parameters
        ----------
        obj : :class:`np.ndarray`
            Rasterized object to be added to the mesh.
        index : :class:`float`
            Refractive index of the object to be added.

        """
        self.mesh *= (-1 * polygone.raster + 1)

        if polygone.Gradient:
            Grad = polygone.Gradient.Evaluate( self.X, self.Y )
            self.mesh += polygone.raster * Grad

        else:
            self.mesh += polygone.raster * polygone.Index



    def CreateMesh(self):
        """ The method create the RI profile mesh according to the user input.

        """

        self.mesh = np.ones(self.Shape)

        self.X, self.Y = np.mgrid[self.Boundaries[0][0]: self.Boundaries[0][1]: complex(self.Shape[0]),
                                  self.Boundaries[1][0]: self.Boundaries[1][1]: complex(self.Shape[1]) ]

        self.coords = np.vstack((self.X.flatten(), self.Y.flatten())).T

        for object in self.Objects:
            self.rasterize_polygone(object)
            self.add_object_to_mesh(object)

        #self.mesh = gaussian_filter(self.mesh, sigma=self.GConv)


    def Plot(self):
        """ The methode plot the rasterized RI profile.

        """
        self.CreateMesh()

        Scene = Scene2D(nCols=1, nRows=1)

        Scene.AddMesh(Row      = 0,
                      Col      = 0,
                      x        = self.X,
                      y        = self.Y,
                      Scalar   = self.mesh,
                      ColorMap = 'coolwarm',
                      xLabel   = r'X-distance [$\mu$m]',
                      yLabel   = r'Y-distance [$\mu$m]',
                      )

        Scene.SetAxes(0, 0, xLimits='Auto', yLimits='Auto', Equal=True)

        Scene.Show()


    def Gradient(self, Plot=False):

        #blurred = gaussian_filter(self.mesh, sigma=0)

        Ygrad, Xgrad = gradientO4( self.mesh.T**2, self.Axes.dx, self.Axes.dy )

        gradient = (Xgrad * self.Axes.XX + Ygrad * self.Axes.YY)

        return gradient.T


class BasedFused():
    def __init__(self, Radius, Fusion, Index, debug, Gradient):
        Mlogger.setLevel(getattr(logging, debug))

        self.Fusion   = Fusion
        self.Radius   = Radius
        self.Index    = Index
        self.Gradient = None

        self.LoadFile()

    def Scale(self, Factor: float):
        self.Object = affinity.scale(self.Object, xfact=Factor, yfact=Factor)

    def LoadFile(self):
        FileName = os.path.join(RootPath, f"Data/Structures/Struc={self.Type}-Fusion={self.Fusion:.1f}")

        with open(FileName, "rb") as poly_file:
            self.C, self.Object, self.hole, self.Scale = pickle.load(poly_file)

    def Plot(self):
        Scene = Scene2D(nCols=1, nRows=1)
        Scene.AddShapely(0, 0, Object=self.Object)

        for n, c in enumerate(self.C):
            Scene.AddShapely(0, 0, Object=c, Text=f"P{n}")

        Scene.SetAxes(0, 0, xLimits='auto', yLimits='auto', Equal=True)
        Scene.Show()


class Circle_(BasedFused):
    def __init__(self, Radius, Index=None, debug='INFO', Gradient=None, Fusion=None, Position=[0,0]):
        self.Type   = 1
        super().__init__(Radius=Radius, Fusion=Fusion, Index=Index, Gradient=Gradient, debug=debug)


    def LoadFile(self):
        FileName = os.path.join(RootPath, f"Data/Structures/Struc={self.Type}")

        with open(FileName, "rb") as poly_file:
            self.C, self.Object, self.hole, self.Scale = pickle.load(poly_file)


class Circle(BasedFused):
    def __init__(self, Position, Radius, Index=None, debug='INFO', Gradient=None):


        assert not all([Index, Gradient]), "Error, you must either define an Index or a Gradient but not both."
        assert any([Index, Gradient]), "Error, you must at least define an Index or a Gradient."
        Mlogger.setLevel(getattr(logging, debug))

        self.Points   = Position
        self.Radius   = Radius
        self.Index    = Index
        self.hole     = None
        self.Gradient = Gradient
        self.Object   = Point(self.Points).buffer(self.Radius)
        self.C        = [Point(Position)]

    def Plot(self):
        PlotObject(Full = [self.Object] + self.C)




class Fused1(BasedFused):
    def __init__(self, Radius, Index=None, debug='INFO', Gradient=None, Fusion=None, Position=[0,0]):
        self.Type   = 1
        super().__init__(Radius=Radius, Fusion=Fusion, Index=Index, Gradient=Gradient, debug=debug)


    def LoadFile(self):
        FileName = os.path.join(RootPath, f"Data/Structures/Struc={self.Type}")

        with open(FileName, "rb") as poly_file:
            self.C, self.Object, self.hole = pickle.load(poly_file)


class Fused4(BasedFused):
    def __init__(self, Radius, Fusion, Index, debug='INFO', Gradient=None):
        self.Type   = 4
        super().__init__(Radius=Radius, Fusion=Fusion, Index=Index, Gradient=Gradient, debug=debug)


class Fused3(BasedFused):
    def __init__(self, Radius, Fusion, Index, debug='INFO', Gradient=None):
        self.Type   = 3
        super().__init__(Radius=Radius, Fusion=Fusion, Index=Index, Gradient=Gradient, debug=debug)


class Fused2(BasedFused):
    def __init__(self, Radius, Fusion, Index, debug='INFO', Gradient=None):
        self.Type   = 2
        super().__init__(Radius=Radius, Fusion=Fusion, Index=Index, Gradient=Gradient, debug=debug)
