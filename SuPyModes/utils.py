import numpy as np
from shapely.geometry.point import Point

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


def NearestPoints(object0, object1):
     return nearest_points(object0.exterior, object1.exterior)


def MakeCircles(Points, Radi):
    C = []
    for point, radius in zip(Points, Radi):
        C.append( Point(point).buffer(radius) )

    return C


def ObjectUnion(object, *args):
    if len(args) == 0:
        args = object.keys()

    objList = [ object[arg] for arg in args ]

    return cascaded_union(objList)


def PlotObject(objects):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    for object in objects:
        ax.add_patch( PolygonPatch(object, alpha=0.5) )

    ax.set_xlim([-100, 100])
    ax.set_ylim([-100, 100])
    ax.set_aspect('equal')
    plt.show()
