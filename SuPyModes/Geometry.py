
""" standard imports """
import sys
import numpy as np
import logging
from numpy           import pi, cos, sin, sqrt, abs, exp, array, ndarray
from matplotlib.path import Path
import matplotlib.pyplot as plt
import matplotlib        as mpl
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools        import combinations
from scipy.optimize   import minimize_scalar
from shapely.geometry import Point, LineString, MultiPolygon, Polygon
from shapely.geometry.collection import GeometryCollection
from shapely.ops import cascaded_union


""" package imports """
from SuPyModes.Directories import RootPath
from SuPyModes.Special     import Intersection
from SuPyModes.utils       import *



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

    def __init__(self, Objects, Xbound, Ybound, Nx, Ny):
        self.Objects    = ToList(Objects)

        self.Indices = sorted( list( set( [1] + [obj.Index for obj in Objects] ) ) )

        self.Boundaries = [Xbound, Ybound]

        self.Shape = [Nx, Ny]

        self.CreateMesh()


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

        poly = cascaded_union(polygone)

        obj_ext = Path(list( polygone.exterior.coords))

        obj_ext = obj_ext.contains_points(self.coords).reshape(self.Shape)

        return obj_ext


    def add_object_to_mesh(self, obj, index):
        """ Method add a rasterized object to the pre-existing profile mesh
        with its provided refractive index.

        Parameters
        ----------
        obj : :class:`np.ndarray`
            Rasterized object to be added to the mesh.
        index : :class:`float`
            Refractive index of the object to be added.

        """

        self.mesh *= (-1 * obj + 1)

        self.mesh += obj * index


    def CreateMesh(self):
        """ The method create the RI profile mesh according to the user input.

        """

        self.mesh = np.ones(self.Shape)

        xv = np.linspace(*self.Boundaries[0], self.Shape[0])
        yv = np.linspace(*self.Boundaries[1], self.Shape[1])

        x, y = np.meshgrid(xv,yv)

        x, y = x.flatten(), y.flatten()

        self.coords = np.vstack((x,y)).T

        for object in self.Objects:
            obj = self.rasterize_polygone(object.Object)
            self.add_object_to_mesh(obj, object.Index)


        self.info = { 'Xbound': self.Boundaries[0],
                      'Ybound': self.Boundaries[1],
                      'Shape':  self.Shape,
                      'Size':   np.size(self.mesh) }


    def Plot(self):
        """ The methode plot the rasterized RI profile.

        """

        fig = plt.figure(figsize=(5,5))

        ax = fig.add_subplot(111)

        ax.set_title('Rasterized RI profil', fontsize=10)
        ax.set_ylabel(r'Y-distance [$\mu$m]')
        ax.set_xlabel(r'X-distance [$\mu$m]')

        x = np.linspace(*self.Boundaries[0], np.shape(self.mesh)[0])
        y = np.linspace(*self.Boundaries[1], np.shape(self.mesh)[1])

        vmin=sorted(self.Indices)[1]/1.1
        vmax=sorted(self.Indices)[-1]

        pcm = ax.pcolormesh(  x,
                              y,
                              self.mesh,
                              cmap    = plt.cm.coolwarm,
                              norm=colors.LogNorm(vmin=vmin, vmax=vmax),
                              shading='auto'
                              )

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=colors.LogNorm(vmin=vmin, vmax=vmax))

        ax.contour(x, y, self.mesh, levels=self.Indices, colors='k')

        sm._A = []
        cbar = plt.colorbar(sm, ax=ax, cax=cax)

        ax.set_aspect('equal')
        plt.show()



class Circle(object):
    def __init__(self, Position, Radi, Index):
        self.Points   = Position
        self.Radi     = Radi
        self.Index    = Index
        self.Init()

    def Init(self):
        self.Object = MakeCircles([self.Points], [self.Radi])[0]
        self.C      = [Point(0,0)]


    def Plot(self):
        PlotObject(Full = [self.Object] + self.C)


class BaseFused():
    def __init__(self, Radius, Fusion, Angle, Theta, Index):
        self.Index   = Index
        self.Radius  = Radius
        self.Fusion  = Fusion
        self.Angle   = Angle
        self.N       = len(self.Angle)
        self.Theta   = Deg2Rad(Theta)
        self.GetFibers()
        self.GetTopology()



    def OptimizeGeometry(self):
        res = minimize_scalar(self.ComputeCost, bounds=self.Bound, method='bounded')
        logging.info(f'Result Rv = {res.x}')
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

        logging.info(f'Topology: {self.Topology} \t Rv: {Rv:.0f} \t Cost: {Cost}')

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
    def __init__(self, Radius, Fusion, Index):
        super().__init__(Radius  = Radius,
                         Fusion  = Fusion,
                         Angle   = [0, 90, 180, 270],
                         Theta   = 45,
                         Index   = Index)

        if self.Topology == 'concave':
            self.Bound = (Radius*1.5, 4000)

        else:
            self.Bound = (0, 4000)

        self.Object = self.OptimizeGeometry()


    def GetTopology(self):
        CenterTriangle = Polygon(self.C)

        Coupler  = ObjectUnion(self.Fibers)

        hole     = CenterTriangle.difference(Coupler)

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
    def __init__(self, Radius, Fusion, Index):
        super().__init__(Radius  = Radius,
                         Fusion  = Fusion,
                         Angle   = [0, 120, 240],
                         Theta   = 30,
                         Index   = Index)

        if self.Topology == 'concave':
            self.Bound = (Radius*1.5, 4000)

        else:
            self.Bound = (0, 4000)

        self.Object = self.OptimizeGeometry()



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
    def __init__(self, Radius, Fusion, Index):
        super().__init__(Radius  = Radius,
                         Fusion  = Fusion,
                         Angle   = [0, 180],
                         Theta   = 90,
                         Index   = Index)

        if self.Topology == 'concave':
            self.Bound = (Radius*1.5, 4000)

        else:
            self.Bound = (0, 4000)

        self.Object = self.OptimizeGeometry()


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
