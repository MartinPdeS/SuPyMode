#-------------------------Importations------------------------------------------
""" standard imports """
import sys, copy, pickle
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
from SuPyModes.utils                 import SortSuperSet
#-------------------------Importations------------------------------------------



class SuPySolver(object):
    """ This object corresponds to the solutioner.
    It solves the eigenvalues problems for a given geometry.

    """

    def __init__(self, Coupler):

        self.Geometry  = Coupler

        self.profile = Coupler.mesh

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

        self.nk = self.profile**2 * self.Axes.Dual.k**2

        self.nk = np.reshape(self.nk, self.info['Size'], order = 'F')


    def laplacian_sparse(self):
        """
        Function that constructs a sparse matrix that applies the 5-point laplacian discretization.

        calls:
            :call1: .initiate_finit_difference_matrix()
        """

        self.Finit_diff = SuPyFinitdifference(self.Axes)

        self.Finit_diff.laplacian_sparse(self.nk, self.Symmetries[1], self.Symmetries[0])


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


        self.Symmetries = [Xsym, Ysym]
        self.debug = debug
        self.Shape = self.info['Shape']

        self.Nstep, self.ITRf, self.Nsol = Nstep, ITRf, Nsol

        ITRList = np.linspace(ITRi, ITRf, Nstep)

        self.Set = SuperSet(IndexProfile = self.profile,
                            NSolutions   = self.Nsol,
                            ITR          = ITRList,
                            Symmetries   = self.Symmetries)

        self.tolerance, self.error = tolerance, error

        self.Axes = SuPyAxes(Meta=metadata)

        self.Solve(ITRList, Nsol)

        return self.Set


    def PreSolution(self):

        if self.iter == 0:
            v0 =  self.profile/np.sum(self.profile)

            max_n = np.max(self.profile)

            beta_square = copy.copy(-(self.Axes.Dual.k*max_n)**2)

        else:

            beta_square = -(self.Set[0].Beta[0]*self.Axes.ITR)**2

            v0 = self.Set[0].Field[-1].ravel()

        return beta_square, v0



    def Solve(self, iteration_list: list, Nsol=1):

        self.iter = 0

        widgets=[Bar('=', '[',  ']'), ' ', Percentage(),  ' ', ETA()]

        bar = ProgressBar(maxval=self.Nstep, widgets=widgets)

        bar.start()

        for n, value in enumerate(iteration_list):
            print(f"{n}/{len(iteration_list)}")

            self.Axes.Scale(value)

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

        betas = np.real(np.sqrt(-values) / self.Axes.ITR)

        vectors = np.real(vectors)

        for solution in range(Nsol):

            index = betas[solution] / self.Axes.Direct.k

            self.Set[solution].Append(Index  = index,
                                      Beta   = betas[solution],
                                      ITR    = self.Axes.ITR,
                                      Field  = vectors[:,solution].reshape(self.Shape),
                                      Axes   = self.Axes)

        if self.debug:
            fig   = plt.figure(figsize=((Nsol)*3,3))
            spec2 = gridspec.GridSpec(ncols=(Nsol), nrows=3, figure=fig)

            for solution in range(Nsol):
                axes = fig.add_subplot(spec2[0:2,solution])
                axes.pcolormesh(vectors[:,solution].reshape(self.Shape).T, shading='auto')
                axes.set_title(f'neff: {betas[solution] / self.Axes.Direct.k:.4f}', fontsize=7)

            plt.show()




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
            self.profile = pickle.load(handle)



    def save_file(self, dir: str):
        """
        arguments:

            :param dir: Directory to save solver file
            :type dir: str.

        calls:
            :call1: XXX
        """

        self.Fields.to_pickle(dir)



# ---
