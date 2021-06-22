import logging
import numpy             as np
import copy              as cp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from descartes           import PolygonPatch
from shapely.geometry    import Point, LineString, MultiPolygon, Polygon
from shapely.geometry.collection import GeometryCollection
from shapely.ops         import nearest_points, cascaded_union
from numpy               import pi, cos, sin, sqrt, abs, exp, array, ndarray
from shapely.affinity    import rotate
from progressbar         import Bar, Percentage, ETA, ProgressBar

from SuPyModes.Config    import *

def RecomposeSymmetries(Input, Axes):

    if Axes.Symmetries[1] == 1:
        Input = np.concatenate((Input[::-1,:],Input),axis=0)

    if Axes.Symmetries[1] == -1:
        Input = np.concatenate((-Input[::-1,:],Input),axis=0)

    if Axes.Symmetries[0] == 1:
        Input = np.concatenate((Input[:,::-1],Input),axis=1)

    if Axes.Symmetries[0] == -1:
        Input = np.concatenate((-Input[:,::-1],Input),axis=1)

    if Axes.Direct.X is not None and Axes.Symmetries[0] != 0:
        Xaxis = np.concatenate( (-Axes.Direct.X[::-1], Axes.Direct.X) )
    else:
        Xaxis = Axes.Direct.X

    if Axes.Direct.Y is not None and Axes.Symmetries[1] != 0:
        Yaxis = np.concatenate( (-Axes.Direct.Y[::-1], Axes.Direct.Y) )
    else:
        Yaxis = Axes.Direct.Y

    return Input, Xaxis, Yaxis


def CheckSymmetries(SuperMode0, SuperMode1):
    Sym0 = SuperMode0.Geometry.Axes.Symmetries
    Sym1 = SuperMode1.Geometry.Axes.Symmetries

    if Sym0[0] == 0 or Sym1[0] == 0: return True

    if Sym0[1] == 0 or Sym1[1] == 0: return True

    if Sym0[0] == - Sym1[0]: return False

    if Sym0[1] == - Sym1[1]: return False

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
            O = mode.Field[n+1] ** mode.Field[n]
            if O < 0.6:
                logging.debug(f'Overlap mismatch at {n} \t Mode:{nm} \t ITR:{itr} \t Overlap:{O}')
                SuperSet.Plot('Fields', iter=n)
                SuperSet.Plot('Fields', iter=n+1)


def SwapProperties(SuperMode0, SuperMode1, N):
    S0, S1 = SuperMode0, SuperMode1

    for p in PROPERTIES:
        getattr(S0, p)[N] = getattr(S1, p)[N]



def SortSuperSet(SuperSet):
    SuperSet1 = cp.deepcopy( SuperSet )

    Overlap = np.zeros(len( SuperSet.Combinations ) )

    for n, itr in Enumerate( SuperSet.Geometry.ITRList[:-1], msg = 'Sorting modes...' ):

        for i, mode0 in enumerate( SuperSet.SuperModes ):

            for j, mode1 in enumerate( SuperSet1.SuperModes ):
                Overlap[j] = mode0.Slice[n].Overlap( mode1.Slice[n+1] )

            if np.max(Overlap) < 0.5:
                logging.debug(n, i,'New mode swapping is occuring...\n', Overlap, '\n\n\n')

            k = np.argmax(Overlap)

            if i != k :
                SuperSet[i].Slice[n+1] = SuperSet1[k].Slice[n+1]

    return SuperSet


def prePlot(func):
    def inner(*args, **kwargs):
        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(111)
        ax.grid()
        kwargs, ybound = func(*args, **kwargs, fig=fig)

        ax.set_ylabel(kwargs['name'] + kwargs['unit'])
        ax.set_xlabel('ITR')

        if kwargs['xlim'] is not None:
            ax.set_xlim(left  = kwargs['xlim'][0],  right = kwargs['xlim'][1])

        if kwargs['ylim'] is not None:
            if kwargs['ylim'][0] is None: ymin = ybound[0]
            else: ymin = kwargs['ylim'][0]

            ymin = max(1e-8, ymin)

            if kwargs['ylim'][1] is None: ymax = ybound[1]
            else: ymax = kwargs['ylim'][1]

            ax.set_ylim( ymin = ymin, ymax = ymax )

        if kwargs['yscale'] == 'log': ax.set_yscale('log')

        if kwargs['xscale'] == 'log': ax.set_xscale('log')

        plt.legend(fontsize=6, title="Modes", fancybox=True)

        return fig

    return inner


def Multipage(filename, figs=None, dpi=200):
    logging.info(f'Saving results into {filename}...')
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()


def animate_plots(base_directory, fname_prefix):
    """
    This function creates a movie of the plots using ffmpeg

    Args:
        base_directory (str): the directory with the plots.
        fname_prefix (str): the filename for the model run

    Returns:
        none but creates mp4 from pngs in base directory

    Author: FJC
    """
    import subprocess

    # animate the pngs using ffmpeg
    system_call = "ffmpeg -framerate 5 -pattern_type glob -i '"+base_directory+"*.png' -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p "+base_directory+fname_prefix+"_movie.mp4"
    subprocess.call(system_call, shell=True)


def Index2NA(nCore, nClad):
    return np.sqrt(nCore**2 - nClad**2)


def NA2nCore(NA, nClad):
    return np.sqrt(NA**2+nClad**2)


def GetWidgetBar(msg):
    return [msg, Bar('=', '[',  ']'), ' ', Percentage(),  ' ', ETA()]


def Enumerate(iterator, msg=''):
    WidgetBar = [msg, Bar('=', '[',  ']'), ' ', Percentage(),  ' ', ETA()]

    bar = ProgressBar(maxval=len(iterator), widgets=WidgetBar)

    bar.start()

    for n, iteration in enumerate(iterator):
        bar.update(n)
        if n == len(iterator)-1:
            bar.finish()
        yield n, iteration


















# -
