#-------------------------Importations------------------------------------------
""" standard imports """
import itertools, pickle, os, progressbar
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import numpy as np


""" package imports """
from SuPyModes.toolbox.SuPyAxes import SuPyAxes
from SuPyModes.toolbox.utils import get_project_root
root = get_project_root()

#-------------------------Importations------------------------------------------



class SuPyProcessing(object):
    """ This object corresponds to the solutioner.
    It solves the eigenvalues problems for a given geometry.

    """

    def __init__(self, Data):

        self.Fields = Data.Fields

        self.setup_axes()


    def compute_gradient(self):
        """
        Compute the index profile gradient and multiply it to the meshgrid to return
        self.gradient.
        """

        Ygradient, Xgradient = self.gradientO4(self.Fields._metadata['profile'].T**2,
                                               self.Axes.real['dx'],
                                               self.Axes.real['dy']
                                               )

        self.gradient = ( Xgradient * self.Axes.real['xmeshes'].T + \
                          Ygradient * self.Axes.real['ymeshes'].T )


    def compute_coupling(self):
        """
        Compute and store coupling, adiabatic and overlap for each couple of modes.
        Debug mode should be actived in this method.
        """

        self.Axes = SuPyAxes(ITR=1, Fields=self.Fields)

        bar = progressbar.ProgressBar(maxval=self.Fields._metadata['Nstep'], \
              widgets=[progressbar.Bar('=',
                                       '[',
                                       ']'),
                                       ' ',
                                       progressbar.Percentage(),
                                       ' ',
                                       progressbar.ETA()])
        bar.start()

        for iter, ITR in enumerate(self.Fields._metadata['Vlist'][:-1]):

            bar.update(iter)

            self.Axes.apply_ITR(ITR)

            self.Axes.apply_scaling(ITR)

            self.compute_gradient()

            couples = list(itertools.combinations_with_replacement(range(self.Fields._metadata['Nsol']),2))

            for couple in couples:
                if iter == 0: continue

                C, A, O = self.Xavier_Daxhelet_coupling(iter, couple[0],  couple[1])

                self.Fields.loc[(iter, couple[0]),('Coupling', couple[1])] = C

                self.Fields.loc[(iter, couple[0]),('Adiabatic', couple[1])] = A

                self.Fields.loc[(iter, couple[0]),('Overlap', couple[1])] = O


        print(self.Fields.loc[1,'Adiabatic'])
        bar.finish()


    def update_degenerate_factor(self, iter, mode0, mode1):
        """
        Update the value of degenerate_factor to take into account the symmetries
        of the system solved.
        """
        self.degenerate_factor = 1

        if self.Fields.loc[(iter - 1, mode0),('std','x_symmetry')] in [1,-1]:
            self.degenerate_factor *= 2

        if self.Fields.loc[(iter - 1, mode0),('std','y_symmetry')] in [1,-1]:
            self.degenerate_factor *= 2


    def verify_symmetries(self, iter, mode0, mode1):
        """
        Verify the symmetries to compute the coupings, if mode0 and mode1 are
        mutually anti-symmetric there is no coupling so it returns False.

        """
        if self.Fields.loc[iter-1, mode0][('std','x_symmetry')] == 0:
            return True

        if self.Fields.loc[iter-1, mode0][('std','y_symmetry')] == 0:
            return True

        if self.Fields.loc[iter-1, mode0][('std','x_symmetry')] == - self.Fields.loc[iter, mode1][('std','x_symmetry')] :
            return False

        elif self.Fields.loc[iter-1, mode0][('std','y_symmetry')] == - self.Fields.loc[iter, mode1][('std','y_symmetry')] :
            return False

        else:
            return True


    def Xavier_Daxhelet_coupling(self, iter, mode0, mode1):
        """
        This method calculate the coupling factor for a certain pair
        of mode.

        arguments:
            :param iter: iteration to compute coupling.
            :type iter: int.

            :param mode0: First mode of pair.
            :type mode0: int.

            :param mode1: Second mode of pair.
            :type mode1: int.

        returns:
            : param Coupling: coupling factor between the two modes
            : type Coupling: float.

            : param Adiabatic: adiabatic factor between the two modes
            : type Adiabatic: float.

            : param Overlap: the sum of overlap between the two modes
            : type Overlap: float.
        """
        if not self.verify_symmetries(iter, mode0, mode1): return 0, np.inf, 0

        self.update_degenerate_factor(iter, mode0, mode1)

        lamb0 = self.Fields.loc[iter-1, mode0][('std','beta')]

        lamb1 = self.Fields.loc[iter, mode1][('std','beta')]

        term1 = -1/2 * complex(0,1)

        term2 = self.Axes.real['k']**2/(np.sqrt(lamb0 * lamb1))

        term3 = np.abs(1/(lamb0 - lamb1))

        Overlap = self.Fields.loc[iter-1, mode0][('std','Field')].reshape(self.Fields._metadata['shape']) * \
                  self.Fields.loc[iter, mode1][('std','Field')].reshape(self.Fields._metadata['shape'])

        new_term = self.gradient * Overlap

        I = np.trapz( [np.trapz(Ix, dx=1) for Ix in new_term], dx=1)

        Coupling = np.abs(self.degenerate_factor * term1 * term2 * term3 * I)

        Adiabatic = np.abs((lamb0 - lamb1)/Coupling)

        return Coupling, Adiabatic, np.sum(Overlap)



    def load_file(self, dir):
        """
        This method load data in pickle (json style) form (.p) and convert
        to dict.

        arguments:
            :param dir: Directory of the .pkl file to load.
            :type dir: str.
        """

        with open(os.path.join(dir), 'rb') as handle:
            self.Fields = pickle.load(handle)


    def setup_axes(self):
        """
        This method assign the datas uploaded from load_json() method.
        """

        self.Axes = SuPyAxes(ITR=1, Fields=self.Fields)


    def save_file(self, dir):
        """ Save file to dir.
        arguments:

            :param dir: Directory to save solver file
            :type dir: str.
        """

        self.Fields.to_pickle(dir)


    def concatenate_figure(self, iter, mode):

        image = self.Fields.loc[(iter,mode)][('std','Field')].reshape(self.Fields._metadata['shape']).astype('float')

        x_axis, y_axis = self.Axes.real['xvector'], self.Axes.real['xvector']

        if self.Fields.loc[(iter,mode)][('std','y_symmetry')] == 1:
            image = np.concatenate((image[::-1,:],image),axis=0)
            y_axis = np.concatenate((-y_axis[::-1],y_axis),axis=0)

        if self.Fields.loc[(iter,mode)][('std','y_symmetry')] == -1:
            image = np.concatenate((-image[::-1,:],image),axis=0)
            y_axis = np.concatenate((-y_axis[::-1],y_axis),axis=0)

        if self.Fields.loc[(iter,mode)][('std','x_symmetry')] == 1:
            image = np.concatenate((image[:,::-1], image),axis=1)
            x_axis = np.concatenate((-x_axis[::-1],x_axis),axis=0)

        if self.Fields.loc[(iter,mode)][('std','x_symmetry')] == -1:
            image = np.concatenate((-image[:,::-1], image),axis=1)
            x_axis = np.concatenate((-x_axis[::-1],x_axis),axis=0)

        return image, x_axis, y_axis


    def debug_mode(self, mode0, mode1=None):

        for iter, ITR in enumerate(self.Fields._metadata['Vlist']):

            image0, _, _ = self.concatenate_figure(iter, mode0)
            image1, _, _ = self.concatenate_figure(iter, mode1)
            fig = plt.figure()
            ax0 = fig.add_subplot(121)
            ax1 = fig.add_subplot(122)

            if mode1:
                plt.title(self.Fields.loc[(iter, 1),('coupling', 0)])
            ax0.imshow(image0)
            ax1.imshow(image1)
            plt.show()


    def plot_solution(self):
        """

        """

        fig, axes = plt.subplots(nrows=1,
                                 ncols=self.Fields._metadata['Nsol'] + 1,
                                 figsize=(3*self.Fields._metadata['Nsol'],3))

        axes[0].pcolormesh(self.Axes.x_vector,
                           self.Axes.y_vector,
                           self.Fields._metadata['profile'].T,
                           norm=colors.PowerNorm(gamma=5),
                           cmap='RdBu')

        axes[0].set_title( 'ITR:{0:.2f}'.format( self.Axes.ITR ) )
        axes[0].set_xlabel(r'X-axis [$\mu$m]')
        axes[0].set_ylabel(r'Y-axis [$\mu$m]')

        for mode, ax in enumerate(axes[1:].flatten()):

            image, x_axis, y_axis = self.concatenate_figure(mode)

            ax.pcolormesh(x_axis,
                          y_axis,
                          image,
                          cmap='RdBu')

            title0 = 'mode:{0}\n'.format(mode)
            title1 = '(index:{0:.5f})\n'.format( self.Field.loc[(iter,mode)][('std','index')] )
            title2 = "x_symmetry:{0}\n".format( self.Field.loc[(iter,mode)][('std','x_symmetry')] )
            title3 = "y_symmetry:{0}\n".format( self.Field.loc[(iter,mode)][('std','y_symmetry')] )
            title4 = "order:{0}".format( self.Fields.loc[(iter,mode)]['order'] )
            ax.set_title(title0+title1+title2+title3+title4, fontsize=10)


        plt.show()


    #4th order accurate gradient function based on 2nd order version from http://projects.scipy.org/scipy/numpy/browser/trunk/numpy/lib/function_base.py
    def gradientO4(self, f, *varargs):
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

# ---
