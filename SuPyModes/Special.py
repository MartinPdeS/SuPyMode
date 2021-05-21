import numpy as np

from SuPyModes.utils import CheckSymmetries, ToList

def FieldOverlap(SuperMode0, SuperMode1, iter):
    return SuperMode0.Field[iter-1] * SuperMode1.Field[iter]


def GeoGradient(Profile, Axes, iter):
    Ygrad, Xgrad = gradientO4(Profile.T**2,
                              Axes[iter].Direct.dx,
                              Axes[iter].Direct.dy )

    return Xgrad * Axes[iter].Direct.XX.T + Ygrad * Axes[iter].Direct.YY.T


def ModeCoupling(SuperMode0, SuperMode1, k, Profile, iter):

        assert SuperMode0.DegenerateFactor == SuperMode1.DegenerateFactor

        if not CheckSymmetries(SuperMode0, SuperMode1): return 0

        beta0 = SuperMode0.Beta[iter-1]

        beta1 = SuperMode1.Beta[iter]

        Delta = beta0 * beta1

        Coupling = -0.5j * k**2 / np.sqrt(Delta)

        Coupling *= np.abs(1/(beta0 - beta1))

        O = FieldOverlap(SuperMode0, SuperMode1, iter)

        G = GeoGradient(Profile, SuperMode0.Axes, iter)

        I = np.trapz( [np.trapz(Ix, dx=1) for Ix in O*G], dx=1)

        Coupling *= I * SuperMode0.DegenerateFactor

        return np.abs(Coupling)


def ModeAdiabatic(SuperMode0, SuperMode1, k, Profile, iter):

    beta0 = SuperMode0.Beta[iter]

    beta1 = SuperMode1.Beta[iter+1]

    Coupling = ModeCoupling(SuperMode0 = SuperMode0,
                            SuperMode1 = SuperMode1,
                            k          = k,
                            Profile    = Profile,
                            iter       = iter)

    Adiabatic = np.abs((beta0 - beta1)/Coupling)

    return Adiabatic


def Overlap(Objects):
    Objects = ToList(Objects)
    object0 = Objects[0]

    for object in Objects:
        object0 = object0.intersection(object)

    return object0


def Intersection(object0, object1):
    return object0.exterior.intersection(object1.exterior)


#4th order accurate gradient function based on 2nd order version from http://projects.scipy.org/scipy/numpy/browser/trunk/numpy/lib/function_base.py
def gradientO4(f, *varargs):
    """Calculate the fourth-order-accurate gradient of an N-dimensional scalar function.
    Uses central differences on the interior and first differences on boundaries
    to give the same shape.
    Inputs:
      f -- An N-dimensional array giving samples of a scalar function
      varargs -- 0, 1, or N scalars giving the sample distances in each direction
    Outputs:
      N arrays of the same shape as f giving the derivative of f with respect
       to each dimension.
    """
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    if n == 0:
        dx = [1.0]*N
    elif n == 1:
        dx = [varargs[0]]*N
    elif n == N:
        dx = list(varargs)
    else:
        raise SyntaxError("invalid number of arguments")

    # use central differences on interior and first differences on endpoints

    outvals = []

    # create slice objects --- initially all are [:, :, ..., :]
    slice0 = [slice(None)]*N
    slice1 = [slice(None)]*N
    slice2 = [slice(None)]*N
    slice3 = [slice(None)]*N
    slice4 = [slice(None)]*N

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'

    for axis in range(N):
        # select out appropriate parts for this dimension
        out = np.zeros(f.shape, f.dtype.char)

        slice0[axis] = slice(2, -2)
        slice1[axis] = slice(None, -4)
        slice2[axis] = slice(1, -3)
        slice3[axis] = slice(3, -1)
        slice4[axis] = slice(4, None)
        # 1D equivalent -- out[2:-2] = (f[:4] - 8*f[1:-3] + 8*f[3:-1] - f[4:])/12.0
        out[tuple(slice0)] = (f[tuple(slice1)] - 8.0*f[tuple(slice2)] + 8.0*f[tuple(slice3)] - f[tuple(slice4)])/12.0

        slice0[axis] = slice(None, 2)
        slice1[axis] = slice(1, 3)
        slice2[axis] = slice(None, 2)
        # 1D equivalent -- out[0:2] = (f[1:3] - f[0:2])
        out[tuple(slice0)] = (f[tuple(slice1)] - f[tuple(slice2)])

        slice0[axis] = slice(-2, None)
        slice1[axis] = slice(-2, None)
        slice2[axis] = slice(-3, -1)
        ## 1D equivalent -- out[-2:] = (f[-2:] - f[-3:-1])
        out[tuple(slice0)] = (f[tuple(slice1)] - f[tuple(slice2)])


        # divide by step size
        outvals.append(out / dx[axis])

        # reset the slice object in this dimension to ":"
        slice0[axis] = slice(None)
        slice1[axis] = slice(None)
        slice2[axis] = slice(None)
        slice3[axis] = slice(None)
        slice4[axis] = slice(None)

    if N == 1:
        return outvals[0]
    else:
        return outvals
