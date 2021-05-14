#-------------------------Importations------------------------------------------
""" standard imports """
import os, json, itertools
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import pandas as pd
import numpy as np


""" package imports """
from SuPyModes.toolbox.utils import get_project_root
root = get_project_root()
#-------------------------Importations------------------------------------------


class Collector(object):


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


    def compute(self, *args, savefig=False):
        """
        This method makes a plot out the Nth solution representing
        the fibers' propagation modes.

        arguments:
            :param *args: input arguments
            :type *args: str.

        calls:
            :call1: .manually_select_fiber()
        """


        self.solution_subplots()

        if savefig:
            plt.savefig('results/coupling_adiabatic.png',format='png')

        plt.ion()

        plt.show()

        input('<Press enter to close()>\n')

        plt.close()


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


        fig0, ax0 = plt.subplots(1)
        fig1, ax1 = plt.subplots(1)
        fig2, ax2 = plt.subplots(1)
        fig3, ax3 = plt.subplots(1, figsize=(1*self.Fields._metadata['Nstep'],3))

        spec2 = gridspec.GridSpec(ncols=self.Fields._metadata['Nsol'], nrows=3, figure=fig2)


        ax0.set_ylim(self.plot_config['index']['y'])
        ax0.set_xlim(self.plot_config['index']['x'])

        ax1.set_ylim(self.plot_config['coupling']['y'])
        ax1.set_xlim(self.plot_config['coupling']['x'])

        ax2.set_ylim(self.plot_config['adiabatic']['y'])
        ax2.set_xlim(self.plot_config['adiabatic']['x'])


        ax2.set_yscale('log')

        ax0.set_ylabel('Effective index')
        ax0.set_xlabel('ITR')

        ax1.set_xlabel('ITR')
        ax1.set_ylabel('Coupling coefficient')

        ax2.set_xlabel('ITR')
        ax2.set_title('Adiabatic criterion')
        ax2.set_ylabel(r'$\frac{1}{\rho} \frac{d\rho}{dz} \,\, [ \mu m^{-1} ]$')

        for j in range(self.Fields._metadata['Nsol']):

            if self.Fields._metadata['naming']:
                name = self.Fields.loc[(self.Fields._metadata['Nstep']-1,j)][('std','name')]

            else:
                name0, name1 = i,j

            ax0.plot(self.Fields.loc[pd.IndexSlice[:, 0],('std','ITR')].values,
                     self.Fields.loc[pd.IndexSlice[:, j], ('std', 'index')],
                     '--',
                     label = '{0}'.format(name))



        combinations = itertools.combinations(list(range(self.Fields._metadata['Nsol'])),2)

        for couple in combinations:

            i, j = couple[0], couple[1]

            if not self.verify_symmetries(0, i, j): continue

            if (i,j) not in [(0,1),(0,3),(0,2),(1,2),(1,3)]: continue

            if self.Fields._metadata['naming']:
                name0 = self.Fields.loc[(self.Fields._metadata['Nstep']-1,i)][('std','name')]
                name1 = self.Fields.loc[(self.Fields._metadata['Nstep']-1,j)][('std','name')]

            else:
                name0, name1 = i,j

            ax1.plot(self.Fields.loc[pd.IndexSlice[:, 0],('std','ITR')].values,
                     self.Fields.loc[pd.IndexSlice[:, i], ('Coupling', j)].values,
                     '--',
                     label = '{0}-{1}'.format(name0,name1))


            ax2.semilogy(self.Fields.loc[pd.IndexSlice[:, 0],('std','ITR')].values,
                         self.Fields.loc[pd.IndexSlice[:, i], ('Adiabatic', j)].values,
                         '--',
                         label = '{0}-{1}'.format(name0,name1))

        ax0.legend(ncol=self.Fields._metadata['Nsol']//2,)
        ax1.legend(ncol=self.Fields._metadata['Nsol']//2,)
        ax2.legend(ncol=self.Fields._metadata['Nsol']//2,)
        ax0.grid()
        ax1.grid()
        ax2.grid()



        for mode in range(self.Fields._metadata['Nsol']):

            image = self.concatenate_figure(self.Fields._metadata['Nstep']-1, mode)

            axes = fig3.add_subplot(spec2[0:2,mode])

            bounds = np.linspace(np.min(image), np.max(image), 200)

            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

            axes.imshow(image, norm=norm, cmap='bwr')

            if isinstance(self.Fields.loc[(self.Fields._metadata['Nstep']-1,mode)][('std','name')], str):
                name = self.Fields.loc[(self.Fields._metadata['Nstep']-1,mode)][('std','name')]
            else:
                name = mode
            axes.set_title('{0}'.format(name))

            axes.axis('off')



        columns = (u'Wavelength [\u03BCm]',
                   'Tolerance',
                   'Error order',
                   'Nx points',
                   'Ny points')


        data = [[self.Fields._metadata['wavelength'],
                 self.Fields._metadata['tolerance'],
                 self.Fields._metadata['error'],
                 self.Fields._metadata['Nx'],
                 self.Fields._metadata['Ny']]]

        rows = ['values']

        ax3.axis('off')
        the_table = ax3.table(cellText=data,
                              rowLabels=rows,
                              bbox = [0, 0, 1, 0.2],
                              rowColours="c",
                              colLabels=columns)




    def concatenate_figure(self, iter, mode):

        image = self.Fields.loc[(self.Fields._metadata['Nstep']-1,mode)][('std','Field')].reshape(self.Fields._metadata['shape']).astype('float')

        if self.Fields.loc[(iter,mode)][('std','y_symmetry')] == 1:
            image = np.concatenate((image[::-1,:],image),axis=0)


        if self.Fields.loc[(iter,mode)][('std','y_symmetry')] == -1:
            image = np.concatenate((-image[::-1,:],image),axis=0)


        if self.Fields.loc[(iter,mode)][('std','x_symmetry')] == 1:
            image = np.concatenate((image[:,::-1], image),axis=1)


        if self.Fields.loc[(iter,mode)][('std','x_symmetry')] == -1:
            image = np.concatenate((-image[:,::-1], image),axis=1)


        return image




# -
