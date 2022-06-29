
""" standard imports """
import os
import numpy as np
import logging
from dataclasses import dataclass

from matplotlib.path             import Path
from itertools                   import combinations
from shapely.geometry            import Point, box
from shapely                     import affinity
from scipy.ndimage.filters       import gaussian_filter
from shapely.geometry            import Point, MultiPolygon, Polygon
from shapely.geometry.collection import GeometryCollection
from scipy.optimize              import minimize_scalar

""" package imports """
from SuPyMode.Tools.Special           import gradientO4
from SuPyMode.Tools.utils             import ToList, Axes
from SuPyMode.Plotting.Plots          import Scene, Axis, Mesh, Contour, ColorBar
from SuPyMode.Tools.utils             import ObjectUnion
from SuPyMode.Tools.ShapelyUtils      import NearestPoints, GetBoundaries, ObjectIntersection, Rotate

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

class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


@dataclass
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

    Clad: object
    Objects: list
    Xbound: list
    Ybound: list
    Nx: int = 100
    Ny: int = 10
    GConv: float = 0
    BackGroundIndex: float = 1.

    def __post_init__(self):
        self.Objects    = ToList(self.Objects)
        self.Boundaries = [self.Xbound, self.Ybound]
        self.Shape      = [self.Nx, self.Ny]

        self.Axes  = Axes( {'wavelength': 1.0,
                            'Xbound'    : self.Xbound,
                            'Ybound'    : self.Ybound,
                            'Nx'        : self.Nx,
                            'Ny'        : self.Ny } )

        self.CreateBackGround()
        self.GetAllIndex()


    def GetFullMesh(self, LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry):

        FullMesh = self._Mesh

        if BottomSymmetry in [1,-1]:
            FullMesh = np.concatenate((FullMesh[::-1, :], FullMesh), axis=1)


        if TopSymmetry in [1, -1]:
            FullMesh = np.concatenate((FullMesh, FullMesh[::-1, :]), axis=1)


        if RightSymmetry in [1, -1]:
            FullMesh = np.concatenate((FullMesh[...], FullMesh[::-1, :]), axis=0)


        if LeftSymmetry in [1, -1]:
            FullMesh = np.concatenate((FullMesh[::-1, :], FullMesh[...]), axis=0)


        return FullMesh


    @property
    def Mesh(self):
        if self._Mesh is None:
            self.CreateMesh()
        return self._Mesh

    @property
    def AllObjects(self):
        return [self.BackGround, self.Clad] + self.Objects

    @property
    def MaxIndex(self):
        ObjectList = [self.Clad] + self.Objects
        return max( [obj.Index for obj in ObjectList] )[0]

    @property
    def MinIndex(self):
        ObjectList = [self.Clad] + self.Objects
        return min( [obj.Index for obj in ObjectList] )[0]


    @property
    def xMax(self):
        return self.Boundaries[0][1]

    @property
    def xMin(self):
        return self.Boundaries[0][0]

    @property
    def yMax(self):
        return self.Boundaries[1][1]

    @property
    def yMin(self):
        return self.Boundaries[1][0]


    def CreateBackGround(self):
        self.BackGround        = Namespace
        xBound = 5*max(abs(self.xMin), abs(self.xMax))
        yBound = 5*max(abs(self.yMin), abs(self.yMax))
        self.BackGround.Object = box(minx=-xBound,
                                     miny=-yBound,
                                     maxx=+xBound,
                                     maxy=+yBound)

        self.BackGround.Index    = self.BackGroundIndex
        self.BackGround.hole     = None
        self.BackGround.Gradient = None


    def GetAllIndex(self,):
        self.AllIndex = []
        for obj in self.AllObjects:
            self.AllIndex.append(float(obj.Index))



    def Rotate(self, Angle):
        for object in self.AllObjects:
            object.Object = affinity.rotate(object.Object, Angle, (0,0))


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

        self._Mesh[np.where(polygone.raster > 0)] = 0

        if polygone.Gradient:
            Grad = polygone.Gradient.Evaluate( self.X, self.Y )
            self._Mesh += polygone.raster * Grad

        else:
            self._Mesh += polygone.raster * polygone.Index


    def CreateMesh(self):
        """ The method create the RI profile mesh according to the user input.

        """

        self._Mesh = np.zeros(self.Shape)

        self.X, self.Y = np.mgrid[self.xMin: self.xMax: complex(self.Shape[0]),
                                  self.yMin: self.yMax: complex(self.Shape[1]) ]

        self.coords = np.vstack((self.X.flatten(), self.Y.flatten())).T

        for object in [self.BackGround, self.Clad] + self.Objects:
            self.rasterize_polygone(object)
            self.add_object_to_mesh(object)

        self._Mesh = gaussian_filter(self._Mesh, sigma=self.GConv)


    def Plot(self):
        """ The methode plot the rasterized RI profile.

        """

        self.CreateMesh()

        Fig = Scene('SuPyMode Figure', UnitSize=(4,4))
        Colorbar = ColorBar(Discreet=True, Position='right')

        ax = Axis(Row              = 0,
                  Col              = 0,
                  xLabel           = r'x [$\mu m$]',
                  yLabel           = r'y [$\mu m$]',
                  Title            = f'Refractive index structure',
                  Legend           = False,
                  Grid             = False,
                  Equal            = True,
                  Colorbar         = Colorbar,
                  xScale           = 'linear',
                  yScale           = 'linear')

        artist = Mesh(X           = self.X,
                      Y           = self.Y,
                      Scalar      = self._Mesh,
                      ColorMap    = 'Blues',
                      )

        ax.AddArtist(artist)

        Fig.AddAxes(ax)

        Fig.Show()




    def _Gradient(self):

        Ygrad, Xgrad = gradientO4( self._Mesh.T**2, self.Axes.dx, self.Axes.dy )

        return Ygrad, Xgrad


    def Gradient(self, Plot=False):

        #blurred = gaussian_filter(self.mesh, sigma=0)

        Ygrad, Xgrad = gradientO4( self._Mesh.T**2, self.Axes.dx, self.Axes.dy )

        gradient = (Xgrad * self.Axes.XX + Ygrad * self.Axes.YY)

        return gradient.T









@dataclass
class BaseFused():
    Radius: float
    Fusion: float
    Angle: float
    Theta: float
    Index: float
    debug: bool
    Gradient: object = None
    
    def __post_init__(self):
        self.hole    = None
        self.N       = len(self.Angle)
        self.Theta   = np.deg2rad(self.Theta)

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
        Scene = Scene2D(nCols=1, nRows=1)
        Scene.AddShapely(0, 0, Object=self.Object)

        for n, c in enumerate(self.C):
            Scene.AddShapely(0, 0, Object=c, Text=f"P{n}")

        Scene.SetAxes(0, 0, xLimits='auto', yLimits='auto', Equal=True)
        Scene.Show()


    def ComputeDifference(self, Added, Removed):
        return abs(Added.area - Removed.Area)


    def GetDistance(self):
        alpha    = (2 - self.N) * np.pi / ( 2 * self.N)

        self.d = ( 1 + np.cos(alpha) ) - np.sqrt(self.N) * np.cos(alpha)

        self.d =  ( self.Radius - ( self.d * self.Radius ) * self.Fusion)

        self.d *= 1 / ( np.cos(alpha) )


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



class Circle(BaseFused):
    def __init__(self, Position, Radius, Index=None, debug='INFO', Gradient=None):


        assert not all([Index, Gradient]), "Error, you must either define an Index or a Gradient but not both."
        assert any([Index, Gradient]), "Error, you must at least define an Index or a Gradient."
        Mlogger.setLevel(getattr(logging, debug))

        self.Points   = Position
        self.Radius   = Radius
        self.Index    = Index
        self.hole     = None
        self.Gradient = Gradient
        self.Object   = Point(Position).buffer(self.Radius)
        self.C        = [Point(Position)]




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

        self.Topology = 'convex' if MaxAdded.area > Removed.Area else 'concave'


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

        s = self.d * np.cos(self.Theta)

        if self.Topology == 'concave':
            H  = -np.sqrt( (Rv - self.Radius)**2 - s**2 )
            sign = 1

        if self.Topology == 'convex':
            H  = -np.sqrt( (self.Radius + Rv)**2 - s**2 )
            sign = -1

        R  = sign * self.d * np.sin(self.Theta) + H

        Points    = Rotate( Object=Point(0, R), Angle=self.Theta * 180/np.pi )

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

        s = self.d * np.cos(self.Theta)

        if self.Topology == 'concave':
            H  = -np.sqrt( (Rv - self.Radius)**2 - s**2 )
            sign = 1

        if self.Topology == 'convex':
            H  = -np.sqrt( (self.Radius + Rv)**2 - s**2 )
            sign = -1

        R  = sign * self.d * np.sin(self.Theta) + H

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
            H  = -np.sqrt( (Rv - self.Radius)**2 - self.d**2 )

        if self.Topology == 'convex':
            H  = -np.sqrt( (self.Radius + Rv)**2 - self.d**2 )

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
