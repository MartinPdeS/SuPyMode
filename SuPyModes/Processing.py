#-------------------------Importations------------------------------------------
""" standard imports """
import itertools, pickle, os, progressbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import pandas as pd
import numpy as np


""" package imports """
from SuPyModes.toolbox.SuPyAxes import SuPyAxes
from SuPyModes.utils            import RecomposeSymmetries, gradientO4


#-------------------------Importations------------------------------------------

class Post(object):
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

        Ygradient, Xgradient = gradientO4(self.Fields._metadata['profile'].T**2,
                                          self.Axes.real['dx'],
                                          self.Axes.real['dy'] )

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

                C, A, O = self.Coupling(iter, couple[0],  couple[1])

                self.Fields.loc[(iter, couple[0]),('Coupling', couple[1])] = C

                self.Fields.loc[(iter, couple[0]),('Adiabatic', couple[1])] = A

                self.Fields.loc[(iter, couple[0]),('Overlap', couple[1])] = O

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


    def Coupling(self, iter, mode0, mode1):
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

        Symmetries = (self.Fields.loc[(iter,mode)][('std','x_symmetry')],
                      self.Fields.loc[(iter,mode)][('std','y_symmetry')])

        return RecomposeSymmetries(image, Symmetries, x_axis, y_axis)


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





class Visual(object):


    def __init__(self, Fields=None):
        if Fields is not None:
            self.Fields = Fields

        self.plot_config = { "index": {"x": [0,1], "y": [1.4,1.46]},
                             "coupling": {"x": None, "y": None},
                             "adiabatic": {"x": [0,1], "y": [1e-4,1e-0]}}


    def coupler_adiabatic(self, length = 1000):
        """
        returns:
            : param dLogRho:
            : type dLogRho: list.

        """

        log_rho = np.log(self.Fields._metadata['Vlist'])

        dLogRho = np.abs(np.convolve(log_rho, [-1,1],'same') / length)

        return dLogRho


    def verify_symmetries(self, iter, mode0, mode1):

        if self.Fields.loc[iter, mode0][('std','x_symmetry')] == 0:
            return True

        if self.Fields.loc[iter, mode0][('std','y_symmetry')] == 0:
            return True

        if self.Fields.loc[iter, mode0][('std','x_symmetry')] == - self.Fields.loc[iter, mode1][('std','x_symmetry')] :
            return False

        elif self.Fields.loc[iter, mode0][('std','y_symmetry')] == - self.Fields.loc[iter, mode1][('std','y_symmetry')] :
            return False

        else:
            return True


    def Plot(self):
        """
        Function that traces the graphs/tables corresponding to the different
        input arguments.

        arguments:
            : param subplot_dic : Dictionary containing information for plotting
            graphs "index", "coupling", "adiabatic", "data_frame".
            : type subplot_dict : dict.

        calls:
            :call1: .save_plots()
            :call2: .plot_results()
        """

        self.PlotIndex()

        self.PlotCoupling()

        self.PlotAdiabatic()

        self.PlotModes()

        plt.show()


    def PlotIndex(self):

        fig, ax = plt.subplots(1)

        Nsol = self.Fields._metadata['Nsol']

        ITR  = self.Fields.loc[pd.IndexSlice[:, 0],('std','ITR')].values

        for mode in range(Nsol):
            index = self.Fields.loc[pd.IndexSlice[:, mode], ('std', 'index')]
            ax.plot(ITR, index, label=f"Mode {mode}")

        ax.legend(fontsize=6)

        ax.grid()


    def PlotCoupling(self):

        fig, ax = plt.subplots(1)

        ax.set_xlabel('ITR')

        ax.set_ylabel('Mode coupling')

        ITR  = self.Fields.loc[pd.IndexSlice[:, 0],('std','ITR')].values

        Nsol = self.Fields._metadata['Nsol']

        combinations = itertools.combinations( np.arange(Nsol), 2 )

        for couple in combinations:
            i, j      = couple[0], couple[1]
            if not self.verify_symmetries(0, i, j): continue
            Coupling  = self.Fields.loc[pd.IndexSlice[:, i], ('Coupling', j)].values
            label     = f"Mode {i} - Mode{j}"
            ax.plot(ITR, Coupling, label=label)

        ax.legend(fontsize=6)

        ax.grid()


    def PlotAdiabatic(self):

        fig, ax = plt.subplots(1)

        ax.set_xlabel('ITR')

        ax.set_ylabel(r'Adiabatic criterion $\frac{1}{\rho} \frac{d\rho}{dz} \,\, [ \mu m^{-1} ]$')

        ITR  = self.Fields.loc[pd.IndexSlice[:, 0],('std','ITR')].values

        Nsol = self.Fields._metadata['Nsol']

        combinations = itertools.combinations( np.arange(Nsol), 2 )

        for couple in combinations:
            i, j      = couple[0], couple[1]
            if not self.verify_symmetries(0, i, j): continue
            Adiabatic = self.Fields.loc[pd.IndexSlice[:, i], ('Adiabatic', j)].values
            label     = f"Mode {i} - Mode{j}"
            ax.plot(ITR, Adiabatic, label=label)

        ax.legend(fontsize=6)

        ax.grid()


    def PlotModes(self):

        Nsol = self.Fields._metadata['Nsol']
        fig = plt.figure(figsize=(Nsol,3))

        spec2 = gridspec.GridSpec(ncols=self.Fields._metadata['Nsol'], nrows=3, figure=fig)

        for mode in range(Nsol):
            image = self.concatenate_figure(self.Fields._metadata['Nstep']-1, mode)
            axes = fig.add_subplot(spec2[0:2,mode])
            axes.imshow(image, cmap='bwr')
            axes.axis('off')
            axes.set_title(f"Mode {mode}")


    def concatenate_figure(self, iter, mode):

        image = self.Fields.loc[(self.Fields._metadata['Nstep']-1,mode)][('std','Field')].reshape(self.Fields._metadata['shape']).astype('float')

        symmetries = (self.Fields.loc[(iter,mode)][('std','x_symmetry')],
                      self.Fields.loc[(iter,mode)][('std','y_symmetry')])

        return RecomposeSymmetries(image, symmetries)[0]


# ---
