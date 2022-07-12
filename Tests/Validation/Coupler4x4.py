
import unittest
from pyface.api        import GUI
import inspect


from SuPyMode.Geometry          import Geometry, Fused4, Circle
from SuPyMode.Solver            import SuPySolver
from SuPyMode.Materials         import FusedSilica
from .utils                     import *

PLOTTIME = 1000 # one seconde

class CouplerTestCase(unittest.TestCase, BaseStepTest):

    @TestFactory('Geometry initialisation')
    def step00(self):
        Wavelength = 1.55e-6
        Index = FusedSilica.GetRI(Wavelength)

        Clad = Fused4(Fusion=0.9, Radius = 62.5, Index = Index)

        Core0 = Circle( Position=Clad.C[0], Radius = 4.1, Index = Index+0.005 )
        Core1 = Circle( Position=Clad.C[1], Radius = 4.1, Index = Index+0.005 )
        Core2 = Circle( Position=Clad.C[2], Radius = 4.1, Index = Index+0.005 )
        Core3 = Circle( Position=Clad.C[3], Radius = 4.1, Index = Index+0.005 )

        self.Geo = Geometry(Clad    = Clad,
                            Objects = [Core0, Core1, Core2, Core3],
                            Xbound  = [-140, 140],
                            Ybound  = [-140, 140],
                            Nx      = 50,
                            Ny      = 50)


    @TestFactory('Geometry plottings')
    def step01(self):
        GUI.invoke_after(PLOTTIME, Close)
        self.Geo.Plot()


    @TestFactory('Solver initialisation')
    def step02(self):
        self.Sol = SuPySolver(Coupler=self.Geo, Tolerance=1e-8, MaxIter = 10000)
        self.Sol.CreateSuperSet(Wavelength=1.55, NStep=300, ITRi=1, ITRf=0.05)


    @TestFactory('Adding modes')
    def step03(self):
        self.Sol.AddModes(Sorting    = 'Field',
                          Symmetries = {'Right': 0, 'Left': 0, 'Top': 0, 'Bottom': 0},
                          nMode      = 6,
                          sMode      = 4 )


    @TestFactory('Getting SuperSet')
    def step04(self):
        self.Set = self.Sol.GetSet()


    @TestFactory('Plottings Fields')
    def step05(self):
        GUI.invoke_after(PLOTTIME, Close)
        self.Set.PlotFields([0, -1])


    @TestFactory('Plottings Coupling')
    def step06(self):
        GUI.invoke_after(PLOTTIME, Close)
        self.Set.PlotCoupling()

    @TestFactory('Plottings Adiabatic')
    def step07(self):
        GUI.invoke_after(PLOTTIME, Close)
        self.Set.PlotAdiabatic()


def suite():
    suite = unittest.TestSuite()
    suite.addTest(CouplerTestCase('test_steps'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner(failfast=True)
    runner.run(suite())
