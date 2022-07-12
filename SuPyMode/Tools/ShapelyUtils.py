from shapely     import affinity
from shapely.ops import nearest_points, unary_union
from SuPyMode.Tools.utils import ToList


def NearestPoints(object0, object1):
     return list( nearest_points(object0.exterior, object1.exterior) )


def GetBoundaries(Objects):
    Objects = ToList(Objects)
    return unary_union(Objects).bounds


def ObjectIntersection(Objects):
    Objects = ToList(Objects)
    object0 = Objects[0]

    for object in Objects:
        object0 = object0.intersection(object)

    return object0


def Rotate(Coord = None, Object=None, Angle=0):

    Angle = ToList(Angle)
    rotated = tuple( affinity.rotate(Object, angle, origin = (0,0) ) for angle in Angle )

    if len(rotated) == 1:
        return rotated[0]
    else:
        return rotated
