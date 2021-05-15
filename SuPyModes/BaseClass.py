
import matplotlib.pyplot as plt

from SuPyModes.Config import *


class SetPlots(object):

    def GetFigure(self, Dict):
        plt.figure(figsize=(8,4))
        plt.grid()
        plt.ylabel(Dict['name'] + Dict['unit'])
        plt.xlabel('ITR')


    def PlotIndex(self):
        I = self.Index
        self.GetFigure(CouplingDict)
        for i in range(self.NSolutions):
            print(i)
            plt.plot(self.ITR, I[i], label=f'{i}')
        plt.legend(fontsize=6, title="Modes", fancybox=True)


    def PlotCoupling(self):
        C = self.Coupling
        self.GetFigure(IndexDict)
        for n, (i,j) in enumerate( self.combinations ):
            plt.plot(self.ITR[:-1], C[n], label=f'{i} - {j}')
        plt.legend(fontsize=6, title="Modes", fancybox=True)


    def PlotAdiabatic(self):
        A = self.Adiabatic
        self.GetFigure(AdiabaticDict)
        for n, (i,j) in enumerate( self.combinations ):
            plt.semilogy(self.ITR[:-1], A[n], label=f'{i} - {j}')
        plt.legend(fontsize=6, title="Modes", fancybox=True)


    def PlotFields(self, iter=None):

        fig   = plt.figure(figsize=(self.NSolutions*3,3))
        spec2 = gridspec.GridSpec(ncols=self.NSolutions, nrows=3, figure=fig)

        for mode in range(self.NSolutions):
            Field, xaxis, yaxis = self[mode].FullField(iter)
            axes = fig.add_subplot(spec2[0:2,mode])
            axes.pcolormesh(xaxis, yaxis, Field, shading='auto')
            axes.set_aspect('equal')
            str   = r"$n_{eff}$"
            title = f"Mode {mode} [{str}: {self[mode].Index[iter]:.4f}]"
            axes.set_title(title, fontsize=8)


    def Plot(self, *args, iter=None):
        self.OrderingModes()

        if 'Index' in args: self.PlotIndex()

        if 'Coupling' in args: self.PlotCoupling()

        if 'Adiabatic' in args: self.PlotAdiabatic()

        if 'Fields' in args: self.PlotFields(iter)

        plt.show()


class SetProperties(object):
    @property
    def Index(self):
        I = []
        for i in range(self.NSolutions):
            I.append( self[i].Index )

        return I


    @property
    def Beta(self):
        B = []
        for i in range(self.NSolutions):
            B.append( self[i].Beta )

        return B


    @property
    def Coupling(self):
        C = []
        for (i,j) in self.combinations:
            C.append( self[i].GetCoupling(self[j]) )

        return C


    @property
    def Adiabatic(self):
        C = []
        for (i,j) in self.combinations:
            C.append( self[i].GetAdiabatic(self[j]) )

        return C
