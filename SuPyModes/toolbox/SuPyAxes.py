
import numpy as np
import copy


class SuPyAxes(object):


    def __init__(self, Meta):

        self.ITR = 1

        self.Nx = Meta['Nx']
        self.Ny = Meta['Ny']

        x_vector = np.linspace(*Meta['Xbound'], Meta['Nx'])
        y_vector = np.linspace(*Meta['Ybound'], Meta['Ny'])

        self.Dual = DirectSpace(xvector      = x_vector,
                                yvector    = y_vector,
                                 wavelength = Meta['wavelength'])

        self.Direct = DirectSpace(xvector    = x_vector,
                                  yvector    = y_vector,
                                  wavelength = Meta['wavelength'])


    def apply_wavelength(self, wavelength):

        self.Dual.wavelength   = wavelength
        self.Direct.wavelength = wavelength
        self.Dual.k            = 2 * np.pi / wavelength
        self.Direct.k          = 2 * np.pi / wavelength


    def Scale(self, Scale):

        self.ITR = Scale

        self.Dual.ApplyITR(Scale)




class DirectSpace():
    def __init__(self, xvector, yvector, wavelength):
        self.ITR = 1
        self.X   = xvector
        self.Y   = yvector
        self.x   = xvector[0], xvector[-1]
        self.y   = yvector[0], yvector[-1]
        self.Nx  = len(xvector)
        self.Ny  = len(yvector)
        self.dx  = np.abs( self.X[0] - self.X[1] )
        self.dy  = np.abs( self.Y[0] - self.Y[1] )
        self.dA  = self.dx * self.dy
        self.YY, self.XX = np.meshgrid(xvector, yvector)
        self.wavelength = wavelength
        self.k          = 2 * np.pi / wavelength

    def Scale(self, factor):
        self.x   *= Scale
        self.y   *= Scale
        self.X   *= Scale
        self.Y   *= Scale
        self.XX  *= Scale
        self.YY  *= Scale
        self.dA  *= Scale**2

    def ApplyITR(self, ITR):
        factor           = ITR / self.ITR
        self.ITR         = ITR
        self.wavelength /= factor
        self.k          *= factor















# -
