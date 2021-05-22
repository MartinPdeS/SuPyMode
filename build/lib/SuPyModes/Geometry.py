#-------------------------Importations------------------------------------------
""" standard imports """
import sys
import numpy as np
from matplotlib.path import Path
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from itertools        import combinations
from scipy.optimize   import minimize_scalar
from shapely.geometry import Point, LineString, MultiPolygon, Polygon
from shapely.geometry.collection import GeometryCollection
pi, cos, sin, sqrt = np.pi, np.cos, np.sin, np.sqrt

from shapely.geometry import Polygon
from shapely.geometry.point import Point
from shapely.ops import cascaded_union
from descartes import PolygonPatch


""" package imports"""
from SuPyModes.Directories import RootPath
from SuPyModes.Special     import Overlap, Intersection
from SuPyModes.utils       import ( NearestPoints,
                                    Deg2Rad,
                                    ObjectIntersection,
                                    Rotate,
                                    GetBoundaries,
                                    MakeCircles,
                                    ObjectUnion,
                                    PlotObject,
                                    ToList        )



import logging
from numpy            import pi, cos, sin, sqrt, abs, exp, array, ndarray

#-------------------------Importations------------------------------------------




class Geometry(object):
    def __init__(self, Objects, Xbound, Ybound, Nx, Ny, Xsym, Ysym):


        #Xsym, Ysym      = Ysym, Xsym  # <-------- something is fucked somewhere else!

        self.Xsym       = Xsym

        self.Ysym       = Ysym

        self.Symmetries = [ Xsym, Ysym ]

        self.Objects    = ToList(Objects)

        self.Boundaries = [Xbound, Ybound]

        self.Shape = [Nx, Ny]

        self.CreateMesh()


    def rasterize_polygone(self, polygone, points, Nx, Ny):

        poly = cascaded_union(polygone)

        obj_ext = Path(list( polygone.exterior.coords))

        obj_ext = obj_ext.contains_points(points).reshape(self.Shape)

        return obj_ext


    def add_object_to_mesh(self, obj, index):

        self.mesh *= (-1 * obj + 1)

        self.mesh += obj * index


    def CreateMesh(self):

        self.mesh = np.ones(self.Shape)

        xv = np.linspace(*self.Boundaries[0], self.Shape[0])
        yv = np.linspace(*self.Boundaries[1], self.Shape[1])

        x, y = np.meshgrid(xv,yv)

        x, y = x.flatten(), y.flatten()

        points = np.vstack((x,y)).T

        for object in self.Objects:
            obj = self.rasterize_polygone(object.Object, points, *self.Shape)
            self.add_object_to_mesh(obj, object.Index)


        self.info = { 'Xbound': self.Boundaries[0],
                      'Ybound': self.Boundaries[1],
                      'Shape':  self.Shape,
                      'Size':   np.size(self.mesh) }


    def Plot(self):

        fig = plt.figure(figsize=(5,5))

        ax = fig.add_subplot(111)

        ax.set_title('Rasterized RI profil', fontsize=10)
        ax.set_ylabel(r'Y-distance [$\mu$m]')
        ax.set_xlabel(r'X-distance [$\mu$m]')

        pcm = ax.pcolormesh(
                            np.linspace(*self.Boundaries[0], np.shape(self.mesh)[0]),
                            np.linspace(*self.Boundaries[1], np.shape(self.mesh)[1]),
                            self.mesh,
                            cmap    = 'PuBu_r',
                            shading = 'auto',
                            norm    = colors.LogNorm(vmin = 1.42, vmax=self.mesh.max()))

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


    def Plot(self):

        fig = plt.figure(figsize=(5,5))

        ax = fig.add_subplot(111)

        ax.grid()

        ax.set_ylabel(r'Distance Y-axis [$\mu m$]')
        ax.set_xlabel(r'Distance X-axis [$\mu m$]')
        ax.set_title(r'Geometrical configuration', fontsize=10)


        for object in self.Object:
            ax.add_patch(PolygonPatch(object, alpha=0.5))


        for n, point in enumerate( self.Points ):
            ax.scatter(point.x, point.y)
            ax.text(point.x, point.y, f'P{n+1}')

        minx, miny, maxx, maxy = cascaded_union(self.Object).bounds

        ax.set_xlim( [ minx, maxx ] )
        ax.set_ylim( [ miny, maxy ] )

        plt.show()


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
        print('Result Rv =',res.x)
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
