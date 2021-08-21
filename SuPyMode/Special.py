import numpy               as np
from numpy                 import zeros, ones, array, sqrt, argsort, absolute, dot
from numpy.linalg          import eig
from scipy.sparse.linalg   import eigs as LA
from scipy.sparse          import diags

from SuPyMode.utils import ToList

def FieldOverlap(SuperMode0, SuperMode1, iter):
    return SuperMode0.Slice[iter-1] * SuperMode1.Slice[iter]


def GeoGradient(Profile, Axes, iter):
    Ygrad, Xgrad = gradientO4(Profile.T**2,
                              Axes.Direct.dx,
                              Axes.Direct.dy )

    return Xgrad * Axes.Direct.XX.T + Ygrad * Axes.Direct.YY.T


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
    dx = list(varargs)

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


def Norm(v):
    return sqrt( dot(v,v) )


def Coefficient(v1, v2):
    return dot(v2, v1) / dot(v1, v1)


def Projection(v1, v2):
    return Coefficient(v1, v2) * v1


def Gram_Schmidt(M):
    Y = copy.copy(M)
    for i in range(1, M.shape[1]):
        for j in range(0,i) :
            temp_vec = M[i] - Projection(Y[j], M[i])

        Y[i] = temp_vec
    return Y


def Normalize(v):
    return v / Norm(v)


def _Gram_Schmidt(A):
    A = normalize(A)

    for i in range(1, A.shape[1]):
        Ai = A[:, i]
        for j in range(0, i):
            Aj = A[:, j]
            t = Ai.dot(Aj)
            Ai = Ai - t * Aj
        A[:, i] = normalize(Ai)

    return A


def Product(v0, v1):
    return dot(v0.conj(), v1)


def Lanczos(A, m):

    m     += 10
    if m > A.shape[1]:
        m = A.shape[1]

    n     = A.shape[1]

    beta  = zeros(m+1);
    alpha = zeros(m)

    V = zeros((n, m))
    v = [zeros(n), Normalize(ones(n))]

    for i in range(0, m):
        w = dot(A, v[1]).squeeze()

        alpha[i] = dot(w, v[1])
        w -= alpha[i] * v[1] - beta[i] * v[0]

        w = OrthogonalizeVector(w, V[:, :i])

        beta[i+1] = Norm(w)

        v[0] = v[1]
        v[1] = Normalize(w)

        V[:, i] = v[0]

    T = diags([beta[1:-1], alpha, beta[1:-1]], [-1, 0, 1]).toarray();

    return V, T


def OrthogonalizeVector(w, V):
    for t in range(V.shape[1]):
      tmpa = dot(w, V[:, t])
      if tmpa == 0.0:
        continue
      w -= tmpa * V[:, t]

    return w


def EigenSolve(A, m):

    V, T = Lanczos(A, m)
    d, v = eig(T)

    sortedIndex = argsort(absolute(d))[::-1]
    d           = d[sortedIndex]
    v           = v[:, sortedIndex]

    s = dot(V, v)

    return d[0:m], s[:, 0:m]
