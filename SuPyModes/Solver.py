#-------------------------Importations------------------------------------------
""" standard imports """
import sys, copy, pickle, progressbar
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import copy
import pandas as pd
import numpy as np
from scipy.sparse.linalg import eigs as LA

""" package imports """
from SuPyModes.toolbox.SuPyAxes import SuPyAxes, _SuPyAxes
from SuPyModes.toolbox.LPModes import LP_names
from SuPyModes.toolbox.SuPyFinitDiff import SuPyFinitdifference
from SuPyModes.SuperMode import SuperMode, SuperSet
from SuPyModes.utils            import RecomposeSymmetries, gradientO4
#-------------------------Importations------------------------------------------



class SuPySolver(object):
    """ This object corresponds to the solutioner.
    It solves the eigenvalues problems for a given geometry.

    """

    def __init__(self, coupler):

        self.profile = coupler.mesh

        self.info = coupler.info

        self.vectors = []

        self.pre_solution = {}


    def generate_dataframe(self, wavelength, Nstep, Nsol, ITRf, tolerance, error):

        index = pd.MultiIndex.from_product([np.arange(self.Nstep),np.arange(self.Nsol)], names=['iteration', 'mode'])


        arrays0 = [['std']*8,
                  ['ITR', 'index','order', 'beta', 'Field', 'x_symmetry', 'y_symmetry', 'name']]

        arrays1 = [['Coupling']*self.Nsol, np.arange(self.Nsol)]

        arrays2 = [['Adiabatic']*self.Nsol, np.arange(self.Nsol)]

        arrays3 = [['Overlap']*self.Nsol, np.arange(self.Nsol)]

        tuples0 = list(zip(*arrays0))
        tuples1 = list(zip(*arrays1))
        tuples2 = list(zip(*arrays2))
        tuples3 = list(zip(*arrays3))

        columns0 = pd.MultiIndex.from_tuples(tuples0)
        columns1 = pd.MultiIndex.from_tuples(tuples1)
        columns2 = pd.MultiIndex.from_tuples(tuples2)
        columns3 = pd.MultiIndex.from_tuples(tuples3)

        Columns = columns0.union(columns1)
        Columns = Columns.union(columns2)
        Columns = Columns.union(columns3)

        def foo(x):
            x[('std','ITR')] = self.Vlist[x.index[0][0]]
            return x

        self.Fields = pd.DataFrame(index=index, columns=Columns)

        self.metadata = {
                    'profile': self.profile,
                    'wavelength': wavelength,
                    'Nstep': Nstep,
                    'Nsol': Nsol,
                    'ITRf': ITRf,
                    'tolerance': tolerance,
                    'error': error,
                    'Vlist': self.Vlist,
                        }

        self.metadata.update(self.info)

        self.Fields._metadata = self.metadata



    def compute_nk(self):
        """
        This method compute the value n**2*k**2 which is one of the termes of
        the eigen value probleme.

        calls:
            :call1: .initiate_finit_difference_matrix()
        """

        self.Axes.update_deltas()

        self.nk = self.profile**2 * self.Axes.sim['k']**2

        self.nk = np.reshape(self.nk, self.info['size'], order = 'F')


    def laplacian_sparse(self):
        """
        Function that constructs a sparse matrix that applies the 5-point laplacian discretization.

        calls:
            :call1: .initiate_finit_difference_matrix()
        """

        self.Finit_diff = SuPyFinitdifference(self.Axes)

        self.Finit_diff.laplacian_sparse(self.nk, self.Xsym, self.Ysym)


    def make_pre_solution(self):
        """
        Generate eigenvector v0 to give to the solver in order to accelerate the
        solving process.

        """

        if self.iter == 0:

            self.pre_solution['v0'] =  self.profile/np.sum(self.profile)

            max_n = np.max(self.profile)

            self.pre_solution['beta_square'] = copy.copy(-(self.Axes.sim['k']*max_n)**2)

        else:

            self.pre_solution['beta_square'] = -(self.Fields.loc[(self.iter-1,0)][('std', 'beta')]*self.Axes.ITR)**2

            self.pre_solution['v0'] = self.vectors[:,0]


    def initiate_finit_difference_matrix(self):
        """
        Generate the finit difference matrix for eigenvalues optimisation
        protocol.

        """

        self.compute_nk()

        self.laplacian_sparse()


    def eigen_sol(self, Nsol: int):
        """
        This method solves the eigenvalues system.

        arguments:
            :param N_solution: Number of solutions of the system we want
            :type N_solution: int.

        """
        print('##########', self.pre_solution['beta_square'])
        self.values, self.vectors = LA(self.Finit_diff.Matrix,
                                            k=Nsol,
                                            M=None,
                                            sigma=self.pre_solution['beta_square'],
                                            which='LM',
                                            v0=self.pre_solution['v0'],
                                            ncv=None,
                                            maxiter=None,
                                            tol=self.tolerance,
                                            return_eigenvectors=True,
                                            Minv=None,
                                            OPinv=None,
                                            OPpart='r'
                                            )

        self.betas = np.real(np.sqrt(-self.values) / self.Axes.ITR)

        self.vectors = np.real(self.vectors)







    def stack_results(self, sol_list: list):

        for solution in sol_list:

            index = self.betas[solution] / self.Axes.real['k']

            self.Fields.loc[(self.iter,solution)][('std','ITR')] = self.Axes.ITR

            self.Fields.loc[(self.iter,solution)][('std','Field')] = self.vectors[:,solution]

            self.Fields.loc[(self.iter,solution)][('std','beta')] = self.betas[solution]

            self.Fields.loc[(self.iter,solution)][('std','index')] = index

            self.Fields.loc[(self.iter,solution)][('std','x_symmetry')] = self.Xsym

            self.Fields.loc[(self.iter,solution)][('std','y_symmetry')] = self.Ysym


    def main_loop(self,
                  wavelength: float,
                  Nstep: int = 2,
                  Nsol: int = 5,
                  Xsym: int = 0,
                  Ysym: int = 0,
                  ITRf: float = 0.1,
                  tolerance: float = 0.0004,
                  error: int = 2,
                  naming: bool = False,
                  debug: bool = False):
        """
        This method compute the solution in a loop that use
        the differents data_set.

        calls:
            :call1: process_solver

        """


        self.Xsym, self.Ysym = Xsym, Ysym

        self.debug, self.naming = debug, naming

        self.Nstep, self.ITRf, self.Nsol = Nstep, ITRf, Nsol

        self.Vlist = np.linspace(1,ITRf, Nstep)

        self.tolerance, self.error = tolerance, error

        self.generate_dataframe(wavelength, Nstep, Nsol, ITRf, tolerance, error)

        self.Axes = SuPyAxes(ITR=1, Fields=self.Fields)

        self.serial_solving(self.Vlist, Nsol)

        self.ordering_data()































    def GetModes(self,
                  wavelength: float,
                  Nstep: int = 2,
                  Nsol: int = 5,
                  Xsym: int = 0,
                  Ysym: int = 0,
                  ITRf: float = 0.1,
                  tolerance: float = 0.0004,
                  error: int = 2,
                  naming: bool = False,
                  debug: bool = False):




        metadata = {'profile': self.profile,
                    'wavelength': wavelength,
                    'Nstep': Nstep,
                    'Nsol': Nsol,
                    'ITRf': ITRf,
                    'tolerance': tolerance,
                    'error': error,
                    'Vlist': [] }

        metadata.update(self.info)

        self.Shape = self.info['shape']

        self.Xsym, self.Ysym = Xsym, Ysym

        self.Nstep, self.ITRf, self.Nsol = Nstep, ITRf, Nsol

        ITRList = np.linspace(1,ITRf, Nstep)

        self.Set = SuperSet(IndexProfile = self.profile,
                            NSolutions   = self.Nsol,
                            ITR          = ITRList)



        self.tolerance, self.error = tolerance, error

        self.Axes = _SuPyAxes(ITR=1, Meta=metadata)

        self._Solve(ITRList, Nsol)

        return self.Set


    def _make_pre_solution(self):

        if self.iter == 0:
            self.pre_solution['v0'] =  self.profile/np.sum(self.profile)

            max_n = np.max(self.profile)

            self.pre_solution['beta_square'] = copy.copy(-(self.Axes.sim['k']*max_n)**2)

        else:
            self.pre_solution['beta_square'] = -(self.Set[0].Beta[0]*self.Axes.ITR)**2

            self.pre_solution['v0'] = self.Set[0].Field[-1].ravel()



    def _Solve(self, iteration_list: list, Nsol=1):

        self.iter = 0

        for _, value in enumerate(iteration_list):

            self.Axes.apply_ITR(value)

            self.initiate_finit_difference_matrix()

            self._make_pre_solution()

            self._eigen_sol(Nsol)

            self.iter += 1


    def _eigen_sol(self, Nsol=1, Tolerance=1e-8, MaxIter=None):

        values, vectors = LA(self.Finit_diff.Matrix,
                             k                   = Nsol,
                             M                   = None,
                             sigma               = self.pre_solution['beta_square'],
                             which               = 'LM',
                             v0                  = self.pre_solution['v0'],
                             ncv                 = None,
                             maxiter             = MaxIter,
                             tol                 = Tolerance,
                             return_eigenvectors = True,
                             Minv                = None,
                             OPinv               = None,
                             OPpart              = 'r'
                             )

        betas = np.real(np.sqrt(-values) / self.Axes.ITR)

        vectors = np.real(vectors)

        for solution in range(Nsol):

            index = betas[solution] / self.Axes.real['k']

            self.Set[solution].Append(Index  = index,
                                        Beta   = betas[solution],
                                        ITR    = self.Axes.ITR,
                                        Field  = vectors[:,solution].reshape(self.Shape),
                                        xSym   = self.Xsym,
                                        ySym   = self.Ysym,
                                        Axes   = self.Axes)




































    def serial_solving(self, iteration_list: list, Nsol=1):
        """
        For loop using the C++ solver to retrieve the eigenvectors and
        eigenvalues

        arguments:
            : param value : XXX
            : type value : float.

        """

        self.iter = 0

        bar = progressbar.ProgressBar(maxval=self.Nstep, \
            widgets=[progressbar.Bar('=',
                                     '[',
                                     ']'),
                                     ' ',
                                     progressbar.Percentage(),
                                     ' ',
                                     progressbar.ETA()])

        bar.start()

        for _, value in enumerate(iteration_list):

            self.Axes.apply_ITR(value)

            self.initiate_finit_difference_matrix()

            self.make_pre_solution()

            self.eigen_sol(Nsol)

            self.stack_results(range(Nsol))

            if self.debug:
                self.plot_solution_debug()

            bar.update(self.iter)

            self.iter += 1

        bar.finish()



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



    def ordering_data(self):

        self.add_order()

        self.Fields.reset_index(level=3, drop=True, inplace=True)

        self.Fields.set_index('order', append=True, inplace=True)

        self.Fields.reset_index(level=1, inplace=True)

        self.Fields.reset_index(level=1, inplace=True)

        self.Fields.index.rename('mode', level=1, inplace=True)


    def add_order(self):

        self.Fields.set_index(('std','x_symmetry'), append=True, inplace=True)
        self.Fields.set_index(('std','y_symmetry'), append=True, inplace=True)

        self.Fields = self.Fields.swaplevel(2,1)
        self.Fields = self.Fields.swaplevel(3,2)

        for iter, row in self.Fields.groupby(level=[0,1,2]):

            temp = -row[('std','beta')].values

            sort = np.argsort(temp).argsort()

            self.Fields.loc[iter,('std','order')] = sort

        final_sort = (-self.Fields.loc[self.Nstep-1][('std','index')]).values.argsort().argsort()

        def add_order(df):
            df['order'] = final_sort
            return df

        def add_name(df):
            x = df.index.get_level_values(1)[0]
            y = df.index.get_level_values(2)[0]
            L = len(df)
            df['std','name'] = LP_names[(x,y)][:L]
            return df

        self.Fields = self.Fields.groupby(level=[0]).apply(add_order)

        if self.naming:
            self.Fields = self.Fields.groupby(level=[0,1,2]).apply(add_name)
            self.metadata.update({'naming': True})
        else: self.metadata.update({'naming': True})

        self.Fields._metadata = self.metadata



    def concatenate_figure(self, iter, mode):

        image = self.Fields.loc[(iter,mode)][('std','Field')].reshape(self.Fields._metadata['shape']).astype('float')

        x_axis, y_axis = self.Axes.real['xvector'], self.Axes.real['xvector']

        Symmetries = (self.Fields.loc[(iter,mode)][('std','x_symmetry')],
                      self.Fields.loc[(iter,mode)][('std','y_symmetry')])

        return RecomposeSymmetries(image, Symmetries, x_axis, y_axis)


    def plot_mode(self, mode, iter, ax=None):

        image, x_axis, y_axis = self.concatenate_figure(mode, iter)

        ax.pcolormesh(x_axis,
                      y_axis,
                      image,
                      cmap='RdBu')

        title0 = 'mode:{0}\n'.format(mode)

        title1 = '(index:{0:.5f})\n'.format( self.Fields.loc[(iter,mode)][('std','index')] )

        title2 = "x_sym.:{0}; y_sym.:{1}\n".format( self.Fields.loc[(iter,mode)][('std','x_symmetry')],
                                                    self.Fields.loc[(iter,mode)][('std','y_symmetry')])

        ax.set_title(title0+title1+title2, fontsize=10)


    def plot_mesh(self, ax):

        ax.pcolormesh(self.Axes.real['xvector'],
                           self.Axes.real['yvector'],
                           self.profile.T,
                           norm=colors.PowerNorm(gamma=5),
                           cmap='RdBu')

        ax.set_title( 'ITR:{0:.2f}'.format( self.Axes.ITR ) )

        ax.set_xlabel(r'X-axis [$\mu$m]')

        ax.set_ylabel(r'Y-axis [$\mu$m]')


    def plot_solution_debug(self):
        """
        arguments:
            :param N_print:
            : type N_print: int.

            :param subplot_dict:
            : type subplot_dict: dict.

        calls:
            :call1: .serial_solving()
        """

        fig, axes = plt.subplots(nrows=1,
                                 ncols=self.Nsol + 1,
                                 figsize=(3*self.Nsol,5)
                                 )

        subset = self.Fields.loc[(self.iter,),]

        subset = subset[subset[('std','beta')].notnull()]

        self.plot_mesh(axes[0])

        axes[0].set_aspect('equal')

        for i, row in subset.iterrows():

            self.plot_mode(i, self.iter, ax=axes[1+i])

            axes[1+i].set_aspect('equal')

            axes[1+i].axis('off')


        plt.show()


# ---
