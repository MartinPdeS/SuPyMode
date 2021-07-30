
""" standard imports """
import sys
import numpy as np
import logging
from numpy                       import pi, cos, sin, sqrt, abs, exp, array, ndarray
from matplotlib.path             import Path
import matplotlib.pyplot         as plt
import matplotlib                as mpl
import matplotlib.colors         as colors
from mpl_toolkits.axes_grid1     import make_axes_locatable
from itertools                   import combinations
from scipy.optimize              import minimize_scalar
from shapely.geometry            import Point, LineString, MultiPolygon, Polygon
from shapely.geometry.collection import GeometryCollection
from shapely.ops                 import cascaded_union
from shapely                     import affinity
from scipy.ndimage.filters import gaussian_filter

""" package imports """
from SuPyModes.Directories       import RootPath
from SuPyModes.Special           import Intersection, gradientO4
from SuPyModes.utils             import *

Mlogger = logging.getLogger(__name__)

class Gradient1:
    def __init__(self, center):
        self.center = center

    def evaluate(self, Xmesh, Ymesh):
        rho = np.sqrt( (Xmesh-self.center[0])**2 + (Ymesh-self.center[1])**2 )+10
        return 1/rho


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

    def __init__(self, Objects, Xbound, Ybound, Nx, Ny, Length=None, debug='INFO'):

        Mlogger.setLevel(getattr(logging, debug))

        self.Objects    = ToList(Objects)

        self.Indices    = sorted( list( set( [1] + [obj.Index for obj in Objects] ) ) )

        self.Boundaries = [Xbound, Ybound]

        self.Shape      = [Nx, Ny]

        self.Length     = Length

        self.Axes  = Axes( {'wavelength': 1.0,
                            'Xbound'    : Xbound,
                            'Ybound'    : Ybound,
                            'Nx'        : Nx,
                            'Ny'        : Ny } )

        self.Symmetries = {'LeftSymmetry'   : 0,
                           'RightSymmetry'  : 0,
                           'TopSymmetry'    : 0,
                           'BottomSymmetry' : 0}



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

        if polygone.Gradient:
            grad = polygone.Gradient(center=(0,0))
            Grad = grad.evaluate(self.X, self.Y)


        self.mesh *= (-1 * polygone.raster + 1)

        self.mesh += polygone.raster * polygone.Index


    def CreateMesh(self):
        """ The method create the RI profile mesh according to the user input.

        """

        self.mesh = np.ones(self.Shape)

        xv = np.linspace(*self.Boundaries[0], self.Shape[0])
        yv = np.linspace(*self.Boundaries[1], self.Shape[1])


        self.X, self.Y = np.meshgrid(yv,xv)

        x, y = self.X.flatten(), self.Y.flatten()

        self.coords = np.vstack((x,y)).T

        for object in self.Objects:
            self.rasterize_polygone(object)
            self.add_object_to_mesh(object)

        self.mesh = gaussian_filter(self.mesh, sigma=0.8)


    def __plot__(self, ax):

        self.CreateMesh()

        Field, xaxis, yaxis = RecomposeSymmetries(self.mesh, self.Symmetries, self.Axes)

        pcm = ax.pcolormesh(  xaxis,
                              yaxis,
                              np.abs(Field),
                              cmap    = plt.cm.coolwarm,
                              shading='auto'
                              )

        ax.set_ylabel(r'Y-distance [$\mu$m]', fontsize=6)

        ax.set_xlabel(r'X-distance [$\mu$m]', fontsize=6)

        divider = make_axes_locatable(ax)

        cax = divider.append_axes("right", size="5%", pad=0.05)

        ax.set_title('Rasterized RI profil', fontsize=10)

        ax.set_aspect('equal')


    def Plot(self):
        """ The methode plot the rasterized RI profile.

        """

        fig = plt.figure(figsize=(5,5))

        ax = fig.add_subplot(111)

        self.__plot__(ax)

        plt.tight_layout()

        plt.show()


    def _Gradient(self):

        Ygrad, Xgrad = gradientO4( self.mesh.T**2, self.Axes.dx, self.Axes.dy )

        return Ygrad, Xgrad


    def Gradient(self, Plot=False):

        #blurred = gaussian_filter(self.mesh, sigma=0)

        Ygrad, Xgrad = gradientO4( self.mesh.T**2, self.Axes.dx, self.Axes.dy )

        gradient = (Xgrad * self.Axes.XX + Ygrad * self.Axes.YY)

        if Plot:
            plt.figure()
            plt.pcolormesh(gradient)
            plt.colorbar()
            plt.show()

        return gradient





class Circle(object):
    def __init__(self, Position, Radi, Index, debug='INFO', Gradient=None):

        Mlogger.setLevel(getattr(logging, debug))

        self.Points   = Position
        self.Radi     = Radi
        self.Index    = Index
        self.hole     = None
        self.Gradient = Gradient
        self.Init()

    def Init(self):
        self.Object = MakeCircles([self.Points], [self.Radi])[0]
        self.C      = [Point(0,0)]


    def Plot(self):
        PlotObject(Full = [self.Object] + self.C)


class BaseFused():
    def __init__(self, Radius, Fusion, Angle, Theta, Index, debug, Gradient=None):

        Mlogger.setLevel(getattr(logging, debug))

        self.Index   = Index
        self.Radius  = Radius
        self.Fusion  = Fusion
        self.Angle   = Angle
        self.hole    = None
        self.N       = len(self.Angle)
        self.Theta   = Deg2Rad(Theta)
        self.Gradient = Gradient
        self.GetFibers()
        self.GetTopology()



    def OptimizeGeometry(self):
        res = minimize_scalar(self.ComputeCost, bounds=self.Bound, method='bounded')
        Mlogger.info(f'Result Rv = {res.x}')
        return self.BuildCoupler(Rv=res.x)


    def GetHole(self):
        CenterPart = Polygon(self.C)

        Coupler  = ObjectUnion(self.Fibers)

        return CenterPart.difference(Coupler)


    def GetFibers(self):
        self.GetDistance()

        self.GetCenters()

        self.Fibers = []

        for i in range(self.N):
            self.Fibers.append( self.C[i].buffer(self.Radius) )

    def Plot(self):
        PlotObject(Full = [self.Object] + self.C)


    def ComputeDifference(self, Added, Removed):
        return abs(Added.area - Removed.Area)


    def GetDistance(self):
        alpha    = (2 - self.N) * pi / ( 2 * self.N)

        self.d = ( 1 + cos(alpha) ) - sqrt(self.N) * cos(alpha)

        self.d =  ( self.Radius - ( self.d * self.Radius ) * self.Fusion)

        self.d *= 1 / ( cos(alpha) )


    def GetMask(self, Circles):
        bound = GetBoundaries(ObjectUnion(Circles))
        n0 = NearestPoints(self.Fibers[2], Circles[0])[0]
        n1 = NearestPoints(self.Fibers[3], Circles[0])[0]

        factor = 4
        P0 = Point(factor*n0.x, factor*n0.y)
        P1 = Point(factor*n1.x, factor*n1.y)
        P2 = Point(0,0)

        mask0  = Polygon([P0, P1, P2])

        mask   = [ Rotate(Object=mask0, Angle=angle) for angle in self.Angle ]

        mask   = ObjectUnion(mask)

        return mask


    def BuildCoupler(self, Rv):
        Circles, Triangles = self.MakeVirtual(Rv=Rv)

        Added = self.GetAdded(Circles, Triangles)

        Removed = self.GetRemoved()

        Coupler = ObjectUnion(self.Fibers + [Added])
        if isinstance(Coupler, GeometryCollection):
            Coupler = MultiPolygon([P for P in Coupler if not isinstance(P, Point)])

        eps = 0.2
        Coupler = Coupler.buffer(eps).buffer(-eps)

        return Coupler


    def ComputeCost(self, Rv):
        Circles, Triangles = self.MakeVirtual(Rv=Rv)

        Added = self.GetAdded(Circles, Triangles)

        Removed = self.GetRemoved()

        Cost = self.ComputeDifference(Added, Removed)

        Mlogger.info(f'Topology: {self.Topology} \t Rv: {Rv:.0f} \t Cost: {Cost}')

        return Cost


    def GetCenters(self):
        P0 = Point(0, self.d)
        Points = Rotate(Object=P0, Angle=self.Angle)
        self.C = [ *Points ]


    def GetAdded(self, Circles, Triangles):
        if self.Topology == 'convex':
            All          = ObjectUnion(Circles+self.Fibers)

            AllTriangles = ObjectUnion(Triangles)

            Added        = AllTriangles.difference(All)

        if self.Topology == 'concave':
            f      = ObjectUnion(self.Fibers)

            b      = ObjectIntersection(Circles)

            edge   = b.difference(f)

            hole   = self.GetHole()

            mask   = self.GetMask(Circles)

            Added  = mask.intersection(edge).difference(hole)

        return Added


    def GetRemoved(self):
        combination = combinations(range(self.N), 2)
        Removed     = []
        for i,j in combination:
            Removed.append( self.Fibers[i].intersection( self.Fibers[j] ) )

        CenterPart = Removed[0]
        for geo in Removed:
            CenterPart = CenterPart.intersection(geo)

        Removed.append(CenterPart)

        area = Removed[0].area * self.N - Removed[-1].area

        Removed = ObjectUnion(Removed)

        Removed.Area = area

        return Removed

    def MakeVirtual(self, Rv):
        Triangles, Points = self.MakeTriangles(Rv=Rv)

        Circles = []
        for point in Points:
            Circles.append(Point(point).buffer(Rv))
            Circles[-1].center = point

        return Circles, Triangles



class Fused4(BaseFused):
    def __init__(self, Radius, Fusion, Index, debug='INFO', Gradient=None):
        super().__init__(Radius  = Radius,
                         Fusion  = Fusion,
                         Angle   = [0, 90, 180, 270],
                         Theta   = 45,
                         Index   = Index,
                         debug   = debug)


        if self.Topology == 'concave':
            self.Bound = (Radius*1.5, 4000)

        else:
            self.Bound = (0, 4000)

        self.Object = self.OptimizeGeometry()
        self.Gradient = Gradient


    def GetTopology(self):
        CenterTriangle = Polygon(self.C)

        Coupler  = ObjectUnion(self.Fibers)

        hole     = CenterTriangle.difference(Coupler)

        self.hole = hole

        Convex   = Coupler.convex_hull.difference(hole)

        MaxAdded = Convex.difference(Coupler)

        Removed  = self.GetRemoved()

        if MaxAdded.area > Removed.Area: self.Topology = 'convex'
        else:                            self.Topology = 'concave'


    def GetAdded(self, Circles, Triangles):
        if self.Topology == 'convex':
            All          = ObjectUnion(Circles+self.Fibers)

            AllTriangles = ObjectUnion(Triangles)

            Added        = AllTriangles.difference(All)

        if self.Topology == 'concave':
            f      = ObjectUnion(self.Fibers)

            b      = ObjectIntersection(Circles)

            edge   = b.difference(f)

            hole   = self.GetHole()

            mask   = self.GetMask(Circles)

            Added  = mask.intersection(edge).difference(hole)

        return Added


    def MakeTriangles(self, Rv):
        self.GetCenters()

        s = self.d * cos(self.Theta)

        if self.Topology == 'concave':
            H  = -sqrt( (Rv - self.Radius)**2 - s**2 )
            sign = 1

        if self.Topology == 'convex':
            H  = -sqrt( (self.Radius + Rv)**2 - s**2 )
            sign = -1

        R  = sign * self.d * sin(self.Theta) + H

        Points    = Rotate( Object=Point(0, R), Angle=self.Theta * 180/pi )

        Points    = [ *Rotate( Object=Points, Angle  = self.Angle ) ]

        Triangle  = Polygon( [ self.C[0], self.C[1], Points[2] ] )

        Triangles = [ *Rotate( Object=Triangle, Angle=self.Angle ) ]

        Triangles = ObjectUnion( Triangles )

        if self.Topology == 'concave':
            Triangles = Rotate( Object=Triangles, Angle=[180] )
            Points    = [ Rotate( Object=pts, Angle=[180] ) for pts in Points]

        return Triangles, Points


class Fused3(BaseFused):
    def __init__(self, Radius, Fusion, Index, debug='INFO', Gradient=None):
        super().__init__(Radius  = Radius,
                         Fusion  = Fusion,
                         Angle   = [0, 120, 240],
                         Theta   = 30,
                         Index   = Index,
                         debug   = debug )

        if self.Topology == 'concave':
            self.Bound = (Radius*1.5, 4000)

        else:
            self.Bound = (0, 4000)

        self.Object = self.OptimizeGeometry()
        self.Gradient = Gradient



    def GetTopology(self):
        hole     = self.GetHole()

        Coupler  = ObjectUnion(self.Fibers)

        Convex   = Coupler.convex_hull.difference(hole)

        MaxAdded = Convex.difference(Coupler)

        Removed  = self.GetRemoved()

        if MaxAdded.area > Removed.Area: self.Topology = 'convex'
        else:                            self.Topology = 'concave'


    def GetMask(self, Circles):
        bound = GetBoundaries(ObjectUnion(Circles))
        n0 = NearestPoints(self.Fibers[0], Circles[0])[0]
        n1 = NearestPoints(self.Fibers[1], Circles[0])[0]

        P0 = Point(2*n0.x, 2*n0.y)
        P1 = Point(2*n1.x, 2*n1.y)
        P2 = Point(0,0)

        mask0  = Polygon([P0, P1, P2])

        mask   = [Rotate(Object=mask0, Angle=angle) for angle in self.Angle]

        mask   = ObjectUnion(mask)

        return mask


    def MakeTriangles(self, Rv):
        self.GetCenters()

        s = self.d * cos(self.Theta)

        if self.Topology == 'concave':
            H  = -sqrt( (Rv - self.Radius)**2 - s**2 )
            sign = 1

        if self.Topology == 'convex':
            H  = -sqrt( (self.Radius + Rv)**2 - s**2 )
            sign = -1

        R  = sign * self.d * sin(self.Theta) + H

        Points    = Rotate( Object=Point(0, R), Angle=240 )

        Points    = [ *Rotate( Object=Points, Angle  = self.Angle ) ]

        Triangle  = Polygon( [ self.C[0], self.C[1], Points[0] ] )

        Triangles = [ *Rotate( Object=Triangle, Angle=self.Angle ) ]

        Triangles = ObjectUnion(Triangles)

        if self.Topology == 'concave':
            Triangles = Rotate( Object=Triangles, Angle=[180] )
            Points    = [ Rotate( Object=pts, Angle=[180] ) for pts in Points]

        return Triangles, Points


class Fused2(BaseFused):
    def __init__(self, Radius, Fusion, Index, debug='INFO', Gradient=None):
        super().__init__(Radius  = Radius,
                         Fusion  = Fusion,
                         Angle   = [0, 180],
                         Theta   = 90,
                         Index   = Index,
                         debug   = debug)

        if self.Topology == 'concave':
            self.Bound = (Radius*1.5, 4000)

        else:
            self.Bound = (0, 4000)

        self.Object = self.OptimizeGeometry()
        self.Gradient = Gradient


    def GetTopology(self):
        Coupler  = ObjectUnion(self.Fibers)

        Convex   = Coupler.convex_hull

        MaxAdded = Convex.difference(Coupler)

        Removed  = self.GetRemoved()

        if MaxAdded.area > Removed.Area: self.Topology = 'convex'
        else:                            self.Topology = 'concave'


    def GetMask(self, Circles):
        bound = GetBoundaries(ObjectUnion(Circles))
        n0 = NearestPoints(self.Fibers[0], Circles[0])[0]
        n1 = NearestPoints(self.Fibers[1], Circles[0])[0]

        factor = 4
        P0 = Point(factor*n0.x, factor*n0.y)
        P1 = Point(factor*n1.x, factor*n1.y)
        P2 = Point(0,0)

        mask0  = Polygon([P0, P1, P2])

        mask   = [ Rotate(Object=mask0, Angle=angle) for angle in self.Angle ]

        mask   = ObjectUnion(mask)

        return mask

    def GetAdded(self, Circles, Triangles):
        if self.Topology == 'convex':
            All          = ObjectUnion(Circles+self.Fibers)

            AllTriangles = ObjectUnion(Triangles)

            Added        = AllTriangles.difference(All)

        if self.Topology == 'concave':
            f      = ObjectUnion(self.Fibers)

            b      = ObjectIntersection(Circles)

            edge   = b.difference(f)

            mask   = self.GetMask(Circles)

            Added  = mask.intersection(edge)


        return Added


    def MakeTriangles(self, Rv):
        self.GetCenters()

        if self.Topology == 'concave':
            H  = -sqrt( (Rv - self.Radius)**2 - self.d**2 )

        if self.Topology == 'convex':
            H  = -sqrt( (self.Radius + Rv)**2 - self.d**2 )

        R  =  abs(H)

        Points    = Rotate( Object=Point(0, R), Angle=90 )

        Points    = [ *Rotate( Object=Points, Angle  = self.Angle ) ]

        Triangle  = Polygon( [ self.C[0], self.C[1], Points[0] ] )

        Triangles = [ *Rotate( Object=Triangle, Angle=self.Angle ) ]

        Triangles = ObjectUnion(Triangles)


        if self.Topology == 'concave':
            Triangles = Rotate( Object=Triangles, Angle=[180] )
            Points    = [ Rotate( Object=pts, Angle=[180] ) for pts in Points]

        return Triangles, Points
