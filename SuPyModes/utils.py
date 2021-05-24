import numpy             as np
import copy              as cp
import logging
import matplotlib.pyplot as plt
from descartes           import PolygonPatch
from shapely.geometry    import Point, LineString, MultiPolygon, Polygon
from shapely.geometry.collection import GeometryCollection
from shapely.ops         import nearest_points, cascaded_union
from numpy               import pi, cos, sin, sqrt, abs, exp, array, ndarray
from shapely.affinity    import rotate

from SuPyModes.Config    import *

def RecomposeSymmetries(Input, Symmetries, Xaxis=None, Yaxis=None):

    if Symmetries[1] == 1:
        Input = np.concatenate((Input[::-1,:],Input),axis=0)

    if Symmetries[1] == -1:
        Input = np.concatenate((-Input[::-1,:],Input),axis=0)

    if Symmetries[0] == 1:
        Input = np.concatenate((Input[:,::-1],Input),axis=1)

    if Symmetries[0] == -1:
        Input = np.concatenate((-Input[:,::-1],Input),axis=1)

    if Xaxis is not None and Symmetries[0] != 0:
        Xaxis = np.concatenate( (-Xaxis[::-1],Xaxis) )

    if Yaxis is not None and Symmetries[1] != 0:
        Yaxis = np.concatenate( (-Yaxis[::-1],Yaxis) )

    return Input, Xaxis, Yaxis


def CheckSymmetries(SuperMode0, SuperMode1):
    if SuperMode0.xSym[0] == 0 or SuperMode1.xSym[0] == 0: return True

    if SuperMode0.ySym[0] == 0 or SuperMode1.ySym[0] == 0: return True

    if SuperMode0.xSym[0] == - SuperMode1.xSym[0]: return False

    if SuperMode0.ySym[0] == - SuperMode1.ySym[0]: return False

    return True


def GetBoundaries(Objects):
    Objects = ToList(Objects)
    return cascaded_union(Objects).bounds


def Overlap(Objects):
    Objects = ToList(Objects)
    object0 = Objects[0]

    for object in Objects:
        object0 = object0.intersection(object)

    return object0


def Deg2Rad(input):
    return input/180*pi


def ObjectUnion(Objects):
    Objects = ToList(Objects)
    return cascaded_union(Objects)


def ObjectIntersection(Objects):
    Objects = ToList(Objects)
    object0 = Objects[0]

    for object in Objects:
        object0 = object0.intersection(object)

    return object0


def NearestPoints(object0, object1):
     return list( nearest_points(object0.exterior, object1.exterior) )


def Rotate(Coord = None, Object=None, Angle=0):

    if Coord:
        x = Coord[0]
        y = Coord[1]
        z = x + 1j*y

        rot = exp(1j*Deg2Rad(Angle))

        z = z * rot
        x = z.real
        y = z.imag

        return (x, y)

    if Object:
        Angle = ToList(Angle)
        rotated = tuple( rotate(Object, angle, origin = (0,0) ) for angle in Angle )
        if len(rotated) == 1: return rotated[0]

        else: return rotated


def PlotGeo(Exterior, Full, ax):
    n = 0
    for object in Exterior:
        if isinstance(object, Point):
            n = PlotPoint(ax, object, n)

        if isinstance(object, LineString):
            ax.plot(*object.coords, linewidth=3)

        if isinstance(object, (Polygon, MultiPolygon) ):
            ax.plot( *object.exterior.xy )

        if isinstance(object, GeometryCollection ):
            for geo in object:
                PlotGeo(Exterior=[geo], Full=[geo], ax=ax)

    for object in Full:
        if isinstance(object, Point):
            n = PlotPoint(ax, object, n)

        if isinstance(object, LineString):
            ax.plot(*object.coords, linewidth=3)

        if isinstance(object, (Polygon, MultiPolygon) ):
            ax.add_patch( PolygonPatch(object, alpha=0.5) )

        if isinstance(object, GeometryCollection ):
            for geo in object:
                PlotGeo(Exterior=[], Full=[geo], ax=ax)


def PlotPoint(ax, point, n):
        ax.scatter(point.x, point.y, linewidth=10)
        ax.text(point.x, point.y, f'P{n}', fontsize=15)
        return n + 1


def PlotObject(Full=[], Exterior=[]):
    Exterior = ToList(Exterior)
    Full     = ToList(Full)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.grid()
    PlotGeo(Exterior, Full, ax)

    factor = 1.3
    minx, miny, maxx, maxy = GetBoundaries(Exterior+Full)
    ax.set_xlim([minx*factor, maxx*factor])
    ax.set_ylim([miny*factor, maxy*factor])
    ax.set_aspect('equal')
    ax.set_ylabel(r'Y-distance [$\mu$m]')
    ax.set_xlabel(r'X-distance [$\mu$m]')
    plt.show()


UlistLike = (list, ndarray, tuple)
def ToList(*args):
    out = []
    for arg in args:
        if not isinstance(arg, UlistLike): out.append( [arg] )
        else: out.append(arg)

    if len(out) == 1: return out[0]
    return out


def MakeCircles(Points, Radi):
    C = []
    for point, radius in zip(Points, Radi):
        C.append( Point(point).buffer(radius) )

    return C


def VerifySorting(SuperSet, Plot=False):
    for nm, mode in enumerate( SuperSet.SuperModes ):
        for n, itr in enumerate( SuperSet.ITR[:-2] ):
            O = GetOverlap(mode.Field[n+1], mode.Field[n])
            if O < 0.6:
                print(f'Error at {n} \t Mode:{nm} \t ITR:{itr} \t Overlap:{O}')
                if Plot:
                    SuperSet.Plot('Fields', iter=n)
                    SuperSet.Plot('Fields', iter=n+1)

def SwapProperties(SuperMode0, SuperMode1, N):
    S0, S1 = SuperMode0, SuperMode1

    for p in PROPERTIES:
        getattr(S0, p)[N] = getattr(S1, p)[N]

def GetOverlap(mode0, mode1):
    return np.abs( np.sum( mode0 * mode1 ) )


def SortSuperSet(SuperSet, Plot=False):
    logging.info('Sorting modes...')
    SuperSet1 = cp.deepcopy(SuperSet)
    Overlap = np.zeros(len( SuperSet.Combinations ) )

    for n, itr in enumerate( SuperSet.ITR[:-1] ):

        for i, mode0 in enumerate( SuperSet.SuperModes ):

            for j, mode1 in enumerate( SuperSet1.SuperModes ):

                Overlap[j] = GetOverlap(mode0.Field[n], mode1.Field[n+1])

            if np.max(Overlap) < 0.5:
                print(n, i,'New mode swapping is occuring...\n', Overlap, '\n\n\n')
                if Plot:
                    SuperSet.Plot('Fields', iter=n)
                    SuperSet1.Plot('Fields', iter=n+1)


            k = np.argmax(Overlap)

            for p in PROPERTIES:
                getattr(SuperSet[i], p)[n+1] = getattr(SuperSet1[k], p)[n+1]

    return SuperSet
