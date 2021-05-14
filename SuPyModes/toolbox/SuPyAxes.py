#-------------------------Importations------------------------------------------
""" standard imports """
import numpy as np
import copy

#-------------------------Importations------------------------------------------

class SuPyAxes(object):


    def __init__(self, ITR, Fields):


        self.lamb = Fields._metadata['wavelength']
        self.ITR = ITR

        self.Nx = Fields._metadata['Nx']
        self.Ny = Fields._metadata['Ny']
        self.k = 2 * np.pi / self.lamb
        x_vector = np.linspace(*Fields._metadata['Xbound'], Fields._metadata['Nx'])
        y_vector = np.linspace(*Fields._metadata['Ybound'], Fields._metadata['Ny'])

        self.real = {'lambda' :  Fields._metadata['wavelength'],
                     'k'      :  2 * np.pi / Fields._metadata['wavelength'],
                     'xbound' :  np.array(Fields._metadata['Xbound']),
                     'ybound' :  np.array(Fields._metadata['Ybound']),
                     'dx'     :  np.abs(Fields._metadata['Xbound'][0]-Fields._metadata['Xbound'][1])/Fields._metadata['Nx'],
                     'dy'     :  np.abs(Fields._metadata['Ybound'][0]-Fields._metadata['Ybound'][1])/Fields._metadata['Ny'],
                     'dA'     :  np.abs(Fields._metadata['Ybound'][0]-Fields._metadata['Ybound'][1])/Fields._metadata['Ny']*np.abs(Fields._metadata['Xbound'][0]-Fields._metadata['Xbound'][1])/Fields._metadata['Nx'],
                     'xvector':  np.linspace(*Fields._metadata['Xbound'], Fields._metadata['Nx']),
                     'yvector':  np.linspace(*Fields._metadata['Ybound'], Fields._metadata['Ny']),
                     'xmeshes' :  np.meshgrid(y_vector, x_vector)[1],
                     'ymeshes' :  np.meshgrid(y_vector, x_vector)[0],
                     }

        self.init = {'lambda' :  Fields._metadata['wavelength'],
                     'k'      :  2 * np.pi / Fields._metadata['wavelength'],
                     'xbound' :  np.array(Fields._metadata['Xbound']),
                     'ybound' :  np.array(Fields._metadata['Ybound']),
                     'dx'     :  np.abs(Fields._metadata['Xbound'][0]-Fields._metadata['Xbound'][1])/Fields._metadata['Nx'],
                     'dy'     :  np.abs(Fields._metadata['Ybound'][0]-Fields._metadata['Ybound'][1])/Fields._metadata['Ny'],
                     'dA'     :  np.abs(Fields._metadata['Ybound'][0]-Fields._metadata['Ybound'][1])/Fields._metadata['Ny']*np.abs(Fields._metadata['Xbound'][0]-Fields._metadata['Xbound'][1])/Fields._metadata['Nx'],
                     'xvector':  np.linspace(*Fields._metadata['Xbound'], Fields._metadata['Nx']),
                     'yvector':  np.linspace(*Fields._metadata['Ybound'], Fields._metadata['Ny']),
                     'xmeshes' :  np.meshgrid(y_vector, x_vector)[1],
                     'ymeshes' :  np.meshgrid(y_vector, x_vector)[0],
                     }

        self.sim = {'lambda' :  Fields._metadata['wavelength'],
                     'k'      :  2 * np.pi / Fields._metadata['wavelength'],
                     'xbound' :  np.array(Fields._metadata['Xbound']),
                     'ybound' :  np.array(Fields._metadata['Ybound']),
                     'dx'     :  np.abs(Fields._metadata['Xbound'][0]-Fields._metadata['Xbound'][1])/Fields._metadata['Nx'],
                     'dy'     :  np.abs(Fields._metadata['Ybound'][0]-Fields._metadata['Ybound'][1])/Fields._metadata['Ny'],
                     'dA'     :  np.abs(Fields._metadata['Ybound'][0]-Fields._metadata['Ybound'][1])/Fields._metadata['Ny']*np.abs(Fields._metadata['Xbound'][0]-Fields._metadata['Xbound'][1])/Fields._metadata['Nx'],
                     'xvector':  np.linspace(*Fields._metadata['Xbound'], Fields._metadata['Nx']),
                     'yvector':  np.linspace(*Fields._metadata['Ybound'], Fields._metadata['Ny']),
                     'xmeshes' :  np.meshgrid(y_vector, x_vector)[1],
                     'ymeshes' :  np.meshgrid(y_vector, x_vector)[0],
                     }




    def apply_wavelength(self, wavelength):
        """
        :param wavelength: value of the wavelength
        :type wavelength: float.

        Function that calculates, from a given wavelength, the corresponding
        k value

        """

        self.sim['lambda'] = wavelength
        self.real['lambda'] = wavelength
        self.sim['k'] = 2 * np.pi / wavelength
        self.real['k'] = 2 * np.pi / wavelength


    def apply_scaling(self, scale=None):

        if not scale:
            scale = self.ITR

        self.real['xbound'] = self.init['xbound'] * scale
        self.real['ybound'] = self.init['ybound'] * scale
        self.real['dx'] = self.init['dx'] * scale
        self.real['dy'] = self.init['dy'] * scale
        self.real['dA'] = self.init['dA'] * scale**2
        self.real['xvector'] = self.init['xvector'] * scale
        self.real['yvector'] = self.init['yvector'] * scale
        self.real['xmeshes'] = self.init['xmeshes'] * scale
        self.real['ymeshes'] = self.init['ymeshes'] * scale


    def apply_ITR(self, ITR):
        """
        :param ITR: value of the ITR
        :type ITR: float.

        """

        self.ITR = ITR

        self.sim['lambda'] = self.real['lambda'] / ITR

        self.sim['k']      = self.real['k'] * ITR


    def update_deltas(self):
        """
        Function that calculates dx, dy and dA the elementary area
        """

        self.sim['dx'] = self.real['dx'] / self.ITR

        self.sim['dy'] = self.real['dy'] / self.ITR

        self.sim['dA'] = self.sim['dx'] * self.sim['dy']





class _SuPyAxes(object):


    def __init__(self, ITR, Meta):


        self.lamb = Meta['wavelength']
        self.ITR = ITR

        self.Nx = Meta['Nx']
        self.Ny = Meta['Ny']
        self.k = 2 * np.pi / self.lamb
        x_vector = np.linspace(*Meta['Xbound'], Meta['Nx'])
        y_vector = np.linspace(*Meta['Ybound'], Meta['Ny'])

        self.real = {'lambda' :  Meta['wavelength'],
                     'k'      :  2 * np.pi / Meta['wavelength'],
                     'xbound' :  np.array(Meta['Xbound']),
                     'ybound' :  np.array(Meta['Ybound']),
                     'dx'     :  np.abs(Meta['Xbound'][0]-Meta['Xbound'][1])/Meta['Nx'],
                     'dy'     :  np.abs(Meta['Ybound'][0]-Meta['Ybound'][1])/Meta['Ny'],
                     'dA'     :  np.abs(Meta['Ybound'][0]-Meta['Ybound'][1])/Meta['Ny']*np.abs(Meta['Xbound'][0]-Meta['Xbound'][1])/Meta['Nx'],
                     'xvector':  np.linspace(*Meta['Xbound'], Meta['Nx']),
                     'yvector':  np.linspace(*Meta['Ybound'], Meta['Ny']),
                     'xmeshes' :  np.meshgrid(y_vector, x_vector)[1],
                     'ymeshes' :  np.meshgrid(y_vector, x_vector)[0],
                     }

        self.init = {'lambda' :  Meta['wavelength'],
                     'k'      :  2 * np.pi / Meta['wavelength'],
                     'xbound' :  np.array(Meta['Xbound']),
                     'ybound' :  np.array(Meta['Ybound']),
                     'dx'     :  np.abs(Meta['Xbound'][0]-Meta['Xbound'][1])/Meta['Nx'],
                     'dy'     :  np.abs(Meta['Ybound'][0]-Meta['Ybound'][1])/Meta['Ny'],
                     'dA'     :  np.abs(Meta['Ybound'][0]-Meta['Ybound'][1])/Meta['Ny']*np.abs(Meta['Xbound'][0]-Meta['Xbound'][1])/Meta['Nx'],
                     'xvector':  np.linspace(*Meta['Xbound'], Meta['Nx']),
                     'yvector':  np.linspace(*Meta['Ybound'], Meta['Ny']),
                     'xmeshes' :  np.meshgrid(y_vector, x_vector)[1],
                     'ymeshes' :  np.meshgrid(y_vector, x_vector)[0],
                     }

        self.sim = {'lambda' :  Meta['wavelength'],
                     'k'      :  2 * np.pi / Meta['wavelength'],
                     'xbound' :  np.array(Meta['Xbound']),
                     'ybound' :  np.array(Meta['Ybound']),
                     'dx'     :  np.abs(Meta['Xbound'][0]-Meta['Xbound'][1])/Meta['Nx'],
                     'dy'     :  np.abs(Meta['Ybound'][0]-Meta['Ybound'][1])/Meta['Ny'],
                     'dA'     :  np.abs(Meta['Ybound'][0]-Meta['Ybound'][1])/Meta['Ny']*np.abs(Meta['Xbound'][0]-Meta['Xbound'][1])/Meta['Nx'],
                     'xvector':  np.linspace(*Meta['Xbound'], Meta['Nx']),
                     'yvector':  np.linspace(*Meta['Ybound'], Meta['Ny']),
                     'xmeshes' :  np.meshgrid(y_vector, x_vector)[1],
                     'ymeshes' :  np.meshgrid(y_vector, x_vector)[0],
                     }




    def apply_wavelength(self, wavelength):
        """
        :param wavelength: value of the wavelength
        :type wavelength: float.

        Function that calculates, from a given wavelength, the corresponding
        k value

        """

        self.sim['lambda'] = wavelength
        self.real['lambda'] = wavelength
        self.sim['k'] = 2 * np.pi / wavelength
        self.real['k'] = 2 * np.pi / wavelength


    def apply_scaling(self, scale=None):

        if not scale:
            scale = self.ITR

        self.real['xbound'] = self.init['xbound'] * scale
        self.real['ybound'] = self.init['ybound'] * scale
        self.real['dx'] = self.init['dx'] * scale
        self.real['dy'] = self.init['dy'] * scale
        self.real['dA'] = self.init['dA'] * scale**2
        self.real['xvector'] = self.init['xvector'] * scale
        self.real['yvector'] = self.init['yvector'] * scale
        self.real['xmeshes'] = self.init['xmeshes'] * scale
        self.real['ymeshes'] = self.init['ymeshes'] * scale


    def apply_ITR(self, ITR):
        """
        :param ITR: value of the ITR
        :type ITR: float.

        """

        self.ITR = ITR

        self.sim['lambda'] = self.real['lambda'] / ITR

        self.sim['k']      = self.real['k'] * ITR


    def update_deltas(self):
        """
        Function that calculates dx, dy and dA the elementary area
        """

        self.sim['dx'] = self.real['dx'] / self.ITR

        self.sim['dy'] = self.real['dy'] / self.ITR

        self.sim['dA'] = self.sim['dx'] * self.sim['dy']
