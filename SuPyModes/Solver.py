
""" standard imports """
import sys, copy, pickle
import logging
import timeit
import copy  as cp
import matplotlib.pyplot             as plt
import matplotlib.gridspec           as gridspec
import numpy                         as np
from progressbar                     import ProgressBar
from scipy.sparse.linalg             import eigsh as LA

""" package imports """
from SuPyModes.toolbox.SuPyAxes      import SuPyAxes
#from SuPyModes.Special               import EigenSolve
from SuPyModes.toolbox.SuPyFinitDiff import SuPyFinitdifference
from SuPyModes.SuperMode             import SuperSet, ModeSlice
from SuPyModes.includes.EigenSolver  import EigenSolving
from SuPyModes.utils                 import SortSuperSet, RecomposeSymmetries, GetWidgetBar, Enumerate
np.set_printoptions(precision=1, linewidth=350)
logging.basicConfig(level=logging.INFO)

Mlogger = logging.getLogger(__name__)


class SuPySolver(object):
    """ This object corresponds to the solutioner.
    It solves the eigenvalues problems for a given geometry.

    """
    def __init__(self, Coupler,  debug='INFO'):
        Mlogger.setLevel(getattr(logging, debug))

        self.Geometry     = Coupler

        self.info         = Coupler.info

        self.vectors      = []

        self.pre_solution = {}

        self._nk          = None

        self._Laplacian   = None


    @property
    def nk(self):
        temp = self.Geometry.mesh**2 * self.Geometry.Axes.Dual.k**2

        self._nk = np.reshape(temp, self.info['Size'], order = 'F')

        return self._nk


    @property
    def Laplacian(self):
        self._Laplacian = SuPyFinitdifference(self.Geometry.Axes)

        self._Laplacian.laplacian_sparse(self.nk,
                                         self.Geometry.Axes.Symmetries[1],
                                         self.Geometry.Axes.Symmetries[0])
        return self._Laplacian


    def GetModes(self,
                  wavelength: float,
                  Nstep:      int   = 2,
                  Nsol:       int   = 5,
                  ITRi:       float = 1.0,
                  ITRf:       float = 0.1,
                  Ysym:       int   = 0,
                  Xsym:       int   = 0,
                  tolerance:  float = 1e-30,
                  error:      int   = 2,
                  naming:     bool  = False,
                  debug:      bool  = False):

        metadata = {'wavelength': wavelength,
                    'Xbound'    : self.info['Xbound'],
                    'Ybound'    : self.info['Ybound'],
                    'Nx'        : self.info['Shape'][0],
                    'Ny'        : self.info['Shape'][1]}


        self.debug = debug

        self.Nstep, self.Nsol = Nstep, Nsol

        self.Geometry.Axes = SuPyAxes(Meta=metadata)

        self.Geometry.ITRList =  np.linspace(ITRi, ITRf, Nstep)

        self.Geometry.Axes.Symmetries = [Xsym, Ysym]

        self.Set = SuperSet(NSolutions = self.Nsol, Geometry = self.Geometry)

        self.tolerance, self.error = tolerance, error

        self.Solve(self.Geometry.ITRList, Nsol)

        return self.Set


    def PreSolution(self):

        if self.iter == 0:
            v0 =  self.Geometry.mesh/np.sum(self.Geometry.mesh)

            max_n = np.max(self.Geometry.mesh)

            beta_square = copy.copy(-(self.Geometry.Axes.Dual.k*max_n)**2)

        else:

            beta_square = -(self.Set[0][0].Beta*self.Geometry.Axes.ITR)**2

            v0 = self.Set[0][-1].ravel()

        return beta_square, v0



    def _Solve(self, iteration_list: list, Nsol=1):
        self.iter = 0


        self.a = EigenSolving(self.Geometry.mesh, 5, 1000, 1e-15)
        self.a.dx = self.Geometry.Axes.Direct.dx
        self.a.dy = self.Geometry.Axes.Direct.dy
        self.a.Lambda = self.Geometry.Axes.Dual.wavelength



        for n, value in Enumerate(iteration_list, msg='Computing super modes: '):

            self.Geometry.Axes.Scale(value)

            self.GetEigenVectors(Nsol)

            self.iter += 1

        self.Set = self.Set.Sort('Fields')



    def Solve(self, iteration_list: list, Nsol=1):
        self.iter = 0

        a = EigenSolving(self.Geometry.mesh, self.Nsol, 100, self.tolerance)

        beta_square, _ = self.PreSolution()

        for n, value in Enumerate(iteration_list, msg='Computing super modes: '):

            self.Geometry.Axes.Scale(value)

            a.dx = self.Geometry.Axes.Direct.dx

            a.dy = self.Geometry.Axes.Direct.dy

            a.Lambda = self.Geometry.Axes.Dual.wavelength

            EigenVectors, EigenValues = a.ComputeEigen(beta_square)

            beta_square, _ = self.PreSolution()

            betas = np.real(np.sqrt(-EigenValues) / self.Geometry.Axes.ITR)

            vectors = np.real(EigenVectors).reshape([Nsol] + self.Geometry.Shape)

            self.AppendSlice(betas, vectors)

            self.iter += 1

        self.Set = self.Set.Sort('Fields')



    def AppendSlice(self, betas, vectors):

        for solution in range(len(betas)):

            index = betas[solution] / self.Geometry.Axes.Direct.k

            Field = vectors[solution,:,:]

            self.Set[solution].Slice.append(ModeSlice( Field, self.Geometry.Axes, Index=index, Beta=betas[solution] ),)



    def GetEigenVectors(self, Nsol=1, MaxIter=None):
        beta_square, v0 = self.PreSolution()

        values, vectors = LA(self.a.GetMatrix(),
                             k                   = Nsol,
                             M                   = None,
                             sigma               = beta_square,
                             which               = 'LM',
                             #v0                  = v0,
                             ncv                 = Nsol*2,
                             maxiter             = None,
                             tol                 = 1e-15,
                             return_eigenvectors = True,
                             Minv                = None,
                             OPinv               = None,
                             #OPpart              = 'r'
                             )

        print(self.a.GetMatrix(),)
        print('sigma',beta_square)
        print('dx',self.Geometry.Axes.Direct.dx)

        betas = np.real(np.sqrt(-values) / self.Geometry.Axes.ITR)

        vectors = np.real(vectors)

        for solution in range(Nsol):
            index = betas[solution] / self.Geometry.Axes.Direct.k

            Field = vectors[:,solution].reshape(self.Geometry.Shape)

            self.Set[solution].Slice.append(ModeSlice( Field, self.Geometry.Axes, Index=index, Beta=betas[solution] ),)




# ---
