import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from SuPyModes.Config import *

from SuPyModes.utils import RecomposeSymmetries, gradientO4

class SuperMode(object):

    def __init__(self, Name, Profile):
        self.Name     = Name
        self.ITR      = []
        self.Index    = []
        self.Field    = []
        self.Beta     = []
        self.xSym     = []
        self.ySym     = []
        self.GeoID    = None
        self.Axes     = []
        self.Profile  = Profile


    def Append(self, **kwargs):
        for property in PROPERTIES:
            getattr(self, property).append( kwargs[property] )


    @property
    def DegenerateFactor(self):
        Factor = 1

        if self.xSym[0] in [1,-1]: Factor *= 2

        if self.ySym[0] in [1,-1]: Factor *= 2

        return Factor


    def Overlap(self, SuperMode, iter):
        Overlap = self.Field[iter-1] * SuperMode.Field[iter]

        return Overlap


    def GeoGradient(self, iter):
        Ygrad, Xgrad = gradientO4(self.Profile.T**2,
                                  self.Axes[iter].real['dx'],
                                  self.Axes[iter].real['dy'] )

        return ( Xgrad * self.Axes[iter].real['xmeshes'].T +
                 Ygrad * self.Axes[iter].real['ymeshes'].T )


    def Coupling(self, SuperMode):
        C = []
        for n, itr in enumerate(self.ITR[:-1]):
            C.append( self._Coupling(SuperMode=SuperMode, iter=n) )

        return C


    def Adiabatic(self, SuperMode):
        A = []
        for n, itr in enumerate(self.ITR[:-1]):
            A.append( self._Adiabatic(SuperMode=SuperMode, iter=n) )

        return A


    def _Coupling(self, SuperMode, iter):

        if not self.CheckSymmetries(SuperMode): return 0

        beta0 = self.Beta[iter-1]

        beta1 = SuperMode.Beta[iter]

        Delta = beta0 * beta1

        Coupling = -0.5j * self.Axes[iter].real['k']**2 / np.sqrt(Delta)

        Coupling *= np.abs(1/(beta0 - beta1))

        term0 = self.Overlap(SuperMode, iter) * self.GeoGradient(iter)

        I = np.trapz( [np.trapz(Ix, dx=1) for Ix in term0], dx=1)

        Coupling *= I * self.DegenerateFactor

        return np.abs(Coupling)


    def _Adiabatic(self, SuperMode, iter):

        beta0 = self.Beta[iter]

        beta1 = SuperMode.Beta[iter+1]

        Coupling = self._Coupling(SuperMode, iter)

        Adiabatic = np.abs((beta0 - beta1)/Coupling)

        return Adiabatic


    def PlotIndex(self):
        plt.figure(figsize=(8,4))
        plt.plot(self.ITR, self.Index)
        plt.xlabel('ITR')
        plt.ylabel('Effective index')
        plt.grid()
        plt.show()

    def PlotBetas(self):
        plt.figure(figsize=(8,4))
        plt.plot(self.ITR, self.Beta)
        plt.xlabel('ITR')
        plt.ylabel(r'$\beta$ factor')
        plt.grid()
        plt.show()

    def PlotCoupling(self, *SuperModes):
        for supermode in SuperModes:
            C = self.Coupling( SuperMode=supermode )

        plt.figure()
        plt.plot(self.ITR[:-1], np.abs(C))
        plt.grid()
        plt.show()


    def PlotAdiabatic(self, *SuperModes):
        for supermode in SuperModes:
            A = self.Adiabatic( SuperMode=supermode )

        plt.figure()
        plt.semilogy(self.ITR[:-1], np.abs(A))
        plt.grid()
        plt.show()


    def CheckSymmetries(self, SuperMode):
        if self.xSym[0] == 0: return True

        if self.ySym[0] == 0: return True

        if self.xSym[0] == - SuperMode.xSym[0]: return False

        if self.ySym[0] == - SuperMode.ySym[0]: return False

        return True


    def PlotFields(self, iter=18):
        plt.figure()

        Field      = self.Field[iter]
        Symmetries = [self.xSym[iter], self.ySym[iter]]

        image, Xaxis, Yaxis = RecomposeSymmetries(Input = Field,
                                                  Symmetries = Symmetries,
                                                  Xaxis = self.Axes[iter].real['xvector'],
                                                  Yaxis = self.Axes[iter].real['yvector'])

        plt.pcolormesh(Xaxis, Yaxis, image.T, shading='auto')

        plt.xlabel(r'Distance x [$\mu$m]')
        plt.ylabel(r'Distance y [$\mu$m]')
        plt.title(self.Name)
        plt.show()



class SuperSet(object):
    def __init__(self, IndexProfile, NSolutions, ITR):
        self.IndexProfile = IndexProfile
        self.NSolutions   = NSolutions
        self.SuperModes   = []
        self.ITR          = ITR
        self.Init()

        self.combinations = tuple(combinations( np.arange(NSolutions), 2 ) )

    def Init(self):
        for solution in range(self.NSolutions):
            supermode = SuperMode(Profile = self.IndexProfile,
                                  Name = f"Mode {solution}")

            self.SuperModes.append(supermode)


    def __getitem__(self, N):
        return self.SuperModes[N]


    @property
    def Index(self):
        I = []
        for i in range(self.NSolutions):
            I.append( self[i].Index )

        return I


    @property
    def Coupling(self):
        C = []
        for (i,j) in self.combinations:
            C.append( self[i].Coupling(self[j]) )

        return C

    @property
    def Adiabatic(self):
        C = []
        for (i,j) in self.combinations:
            C.append( self[i].Adiabatic(self[j]) )

        return C


    def PlotIndex(self):
        I = self.Index

        plt.figure(figsize=(8,4))

        for i in range(self.NSolutions):
            plt.plot(self.ITR, I[i], label=f'{i}')

        plt.grid()
        plt.legend(fontsize=6)
        plt.xlabel('ITR')
        plt.ylabel('Index')


    def PlotCoupling(self):
        C = self.Coupling

        plt.figure(figsize=(8,4))

        for n, (i,j) in enumerate( self.combinations ):
            plt.plot(self.ITR[:-1], C[n], label=f'{i} - {j}')

        plt.grid()
        plt.legend()
        plt.xlabel('ITR')
        plt.ylabel('Mode coupling')


    def PlotAdiabatic(self):
        A = self.Adiabatic

        plt.figure(figsize=(8,4))

        for n, (i,j) in enumerate( self.combinations ):
            plt.semilogy(self.ITR[:-1], A[n], label=f'{i} - {j}')

        plt.grid()
        plt.legend(fontsize=6)
        plt.xlabel('ITR')

        plt.ylabel(AdiabaticDict['name'] + AdiabaticDict['unit'])



    def Plot(self, *args):
        self.OrderingModes()
        if 'Index' in args:
            self.PlotIndex()

        if 'Coupling' in args:
            self.PlotCoupling()

        if 'Adiabatic' in args:
            self.PlotAdiabatic()

        plt.show()


    def SwapProperties(self, SuperMode0, SuperMode1, N):
        S0 = SuperMode0
        S1 = SuperMode1
        for p in PROPERTIES:
            getattr(S0, p)[N], getattr(S1, p)[N] = getattr(S1, p)[N], getattr(S0, p)[N]


    def OrderingModes(self):
        for (i,j) in self.combinations:
            for n, ITR in enumerate(self.ITR):
                if i < j  and self[i].Index[n] > self[j].Index[n]:
                    self.SwapProperties(self[i], self[j], n)



















# -
