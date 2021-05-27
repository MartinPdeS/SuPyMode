#-------------------------Importations------------------------------------------
""" standard imports """
import sys, copy, pickle
import logging
from progressbar           import Bar, Percentage, ETA, ProgressBar
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
import numpy               as np
from scipy.sparse.linalg   import eigs as LA

""" package imports """
from SuPyModes.toolbox.SuPyAxes      import SuPyAxes
from SuPyModes.toolbox.LPModes       import LP_names
from SuPyModes.toolbox.SuPyFinitDiff import SuPyFinitdifference
from SuPyModes.SuperMode             import SuperSet
from SuPyModes.utils                 import SortSuperSet, RecomposeSymmetries
#-------------------------Importations------------------------------------------

import logging
logging.basicConfig(level=logging.INFO)

class SuPySolver(object):
    """ This object corresponds to the solutioner.
    It solves the eigenvalues problems for a given geometry.

    """
    def __init__(self, Coupler):

        self.Geometry  = Coupler

        self.info = Coupler.info

        self.vectors = []

        self.pre_solution = {}



    def compute_nk(self):
        """
        This method compute the value n**2*k**2 which is one of the termes of
        the eigen value probleme.

        calls:
            :call1: .initiate_finit_difference_matrix()
        """

        self.nk = self.Geometry.mesh**2 * self.Geometry.Axes.Dual.k**2

        self.nk = np.reshape(self.nk, self.info['Size'], order = 'F')


    def laplacian_sparse(self):
        """
        Function that constructs a sparse matrix that applies the 5-point laplacian discretization.

        calls:
            :call1: .initiate_finit_difference_matrix()
        """

        self.Finit_diff = SuPyFinitdifference(self.Geometry.Axes)

        self.Finit_diff.laplacian_sparse(self.nk,
                                         self.Geometry.Axes.Symmetries[1],
                                         self.Geometry.Axes.Symmetries[0])


    def initiate_finit_difference_matrix(self):
        """ Generate the finit difference matrix for eigenvalues optimisation
        protocol.

        """
        self.compute_nk()

        self.laplacian_sparse()


    def GetModes(self,
                  wavelength: float,
                  Nstep:      int   = 2,
                  Nsol:       int   = 5,
                  ITRi:       float = 1.0,
                  ITRf:       float = 0.1,
                  Ysym:       int   = 0,
                  Xsym:       int   = 0,
                  tolerance:  float = 1e-16,
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

            beta_square = -(self.Set[0].Beta[0]*self.Geometry.Axes.ITR)**2

            v0 = self.Set[0].Field[-1].ravel()

        return beta_square, v0



    def Solve(self, iteration_list: list, Nsol=1):

        self.iter = 0

        widgets=[Bar('=', '[',  ']'), ' ', Percentage(),  ' ', ETA()]

        bar = ProgressBar(maxval=self.Nstep, widgets=widgets)

        bar.start()

        for n, value in enumerate(iteration_list):
            logging.info(f"{n}/{len(iteration_list)}")

            self.Geometry.Axes.Scale(value)

            self.initiate_finit_difference_matrix()

            self.GetEigenVectors(Nsol)

            self.iter += 1

            bar.update(self.iter)

        self.Set = SortSuperSet(self.Set)




    def GetEigenVectors(self, Nsol=1, MaxIter=None):
        beta_square, v0 = self.PreSolution()

        values, vectors = LA(self.Finit_diff.Matrix,
                             k                   = Nsol,
                             M                   = None,
                             sigma               = beta_square,
                             which               = 'LR',
                             v0                  = v0,
                             ncv                 = None,
                             maxiter             = MaxIter,
                             tol                 = 1e-30,
                             return_eigenvectors = True,
                             Minv                = None,
                             OPinv               = None,
                             OPpart              = 'r' )

        betas = np.real(np.sqrt(-values) / self.Geometry.Axes.ITR)

        vectors = np.real(vectors)

        for solution in range(Nsol):
            index = betas[solution] / self.Geometry.Axes.Direct.k

            Field = vectors[:,solution].reshape(self.Geometry.Shape)

            self.Set[solution].Append(Index  = index,
                                      Beta   = betas[solution],
                                      ITR    = self.Geometry.Axes.ITR,
                                      Field  = ModeArray( Field, self.Geometry.Axes ),
                                      Axes   = self.Geometry.Axes)




    def load_file(self, dir: str):
        """
        This method load data in pickle (json style) form (.p) and convert
        to dict.

        arguments:
            :param dir: Directory of the .p file to load.
            :type dir: str.

        calls:
        :call1:
        """

        with open(dir, 'rb') as handle:
            self.Geometry.mesh = pickle.load(handle)



    def save_file(self, dir: str):
        """
        arguments:

            :param dir: Directory to save solver file
            :type dir: str.

        calls:
            :call1: XXX
        """

        self.Fields.to_pickle(dir)




class ModeArray(np.ndarray):

    def __new__(cls, input_array, Axes=None):
        self = input_array.view(ModeArray)
        self.Axes = Axes
        return self

    def __array_finalize__(self, viewed):
        pass


    def __pow__(self, other):
        assert isinstance(other, ModeArray), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    def __plot__(self, ax, title=None):
        Field, xaxis, yaxis = RecomposeSymmetries(self, self.Axes)

        ax.pcolormesh(yaxis, xaxis, Field, shading='auto')

        ax.set_ylabel(r'Y-distance [$\mu$m]', fontsize=6)

        ax.set_xlabel(r'X-distance [$\mu$m]', fontsize=6)

        ax.set_aspect('equal')
        if title:
            ax.set_title(title, fontsize=8)




# ---
