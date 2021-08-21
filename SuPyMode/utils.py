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

from SuPyMode.Config    import *
import numpy as np




class Axes(object):
    def __init__(self, Meta):
        self._k = None

        self.ITR = 1

        self.Nx = Meta['Nx']
        self.Ny = Meta['Ny']

        self.X = np.linspace(*Meta['Xbound'], Meta['Nx'])
        self.Y = np.linspace(*Meta['Ybound'], Meta['Ny'])

        self.XX, self.YY = np.meshgrid(self.X, self.Y)

        self.rho = np.sqrt(self.XX**2 + self.YY**2)

        self.dx  = np.abs( self.X[0] - self.X[1] )
        self.dy  = np.abs( self.Y[0] - self.Y[1] )

        self._wavelength  = Meta['wavelength']
        self._k           = 2 * np.pi / self._wavelength

        self.dA  = self.dx * self.dy

        self.LeftSymmetry   = None
        self.RightSymmetry  = None
        self.TopSymmetry    = None
        self.BottomSymmetry = None

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        self._wavelength = value
        self._k          = 2*np.pi / value

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, value):
        self._k          = value
        self._wavelength = 2*np.pi / value


def RecomposeSymmetries(Input, Symmetries, Axes):

    Xaxis = Axes.Y
    Yaxis = Axes.X

    if Symmetries['BottomSymmetry'] == 1:
        Input = np.concatenate((Input[:,::-1], Input), axis=1)
        Xaxis = np.concatenate( [Xaxis + Xaxis[0] - Xaxis[-1], Xaxis] )


    elif Symmetries['BottomSymmetry'] == -1:
        Input = np.concatenate((-Input[:,::-1], Input), axis=1)
        Xaxis = np.concatenate( [Xaxis + Xaxis[0] - Xaxis[-1], Xaxis] )


    if Symmetries['TopSymmetry'] == 1:
        Input = np.concatenate((Input, Input[:,::-1]), axis=1)
        Xaxis = np.concatenate( [Xaxis, Xaxis - Xaxis[0] + Xaxis[-1]] )

    elif Symmetries['TopSymmetry'] == -1:
        Input = np.concatenate((Input, -Input[:,::-1]), axis=1)
        Xaxis = np.concatenate( [Xaxis, Xaxis - Xaxis[0] + Xaxis[-1]] )


    if Symmetries['RightSymmetry'] == 1:
        Input = np.concatenate((Input, Input[::-1,:]), axis=0)
        Yaxis = np.concatenate( [Yaxis, Yaxis - Yaxis[0] + Yaxis[-1]] )

    elif Symmetries['RightSymmetry'] == -1:
        Input = np.concatenate((Input, -Input[::-1,:]), axis=0)
        Yaxis = np.concatenate( [Yaxis, Yaxis - Yaxis[0] + Yaxis[-1]] )


    if Symmetries['LeftSymmetry'] == 1:
        Input = np.concatenate((Input[::-1,:], Input), axis=0)
        Yaxis = np.concatenate( [Yaxis + Yaxis[0] - Yaxis[-1], Yaxis] )

    elif Symmetries['LeftSymmetry'] == -1:
        Input = np.concatenate((-Input[::-1,:], Input), axis=0)
        Yaxis = np.concatenate( [Yaxis + Yaxis[0] - Yaxis[-1], Yaxis] )

    return Input, Xaxis, Yaxis


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
