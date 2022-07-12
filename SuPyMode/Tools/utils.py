import logging
import numpy             as np
import copy              as cp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from descartes           import PolygonPatch
from shapely.geometry    import Point, LineString, MultiPolygon, Polygon
from shapely.geometry.collection import GeometryCollection
from shapely.ops         import nearest_points, unary_union
from numpy               import pi, cos, sin, sqrt, abs, exp, array, ndarray
from shapely.affinity    import rotate


from SuPyMode.Tools.Config    import *
import numpy as np




class Axes(object):
    def __init__(self, Meta):

        self.ITR = 1

        self.Nx = Meta['Nx']
        self.Ny = Meta['Ny']

        self.X = np.linspace(*Meta['Xbound'], Meta['Nx'])
        self.Y = np.linspace(*Meta['Ybound'], Meta['Ny'])

        self.XX, self.YY = np.meshgrid(self.X, self.Y)

        self.rho = np.sqrt(self.XX**2 + self.YY**2)

        self.dx  = np.abs( self.X[0] - self.X[1] )
        self.dy  = np.abs( self.Y[0] - self.Y[1] )

        self.dA  = self.dx * self.dy




def GetBoundaries(Objects):
    Objects = ToList(Objects)
    return unary_union(Objects).bounds


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
    return unary_union(Objects)


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


def ToList(*args):
    out = []
    for arg in args:
        if not isinstance(arg, (list, ndarray, tuple)): out.append( [arg] )
        else: out.append(arg)

    if len(out) == 1: return out[0]
    return out


def MakeCircles(Points, Radi):
    C = []
    for point, radius in zip(Points, Radi):
        C.append( Point(point).buffer(radius) )

    return C




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
















# -
