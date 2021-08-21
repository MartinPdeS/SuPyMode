from SuPyModes.Geometry          import Geometry, Circle, Fused2, Gradient
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
from SuPyModes.fibers            import *
import matplotlib.pyplot         as plt

"""
FIGURE 2.5 SBB_____________________________________________________

"""

A = Fiber_DCF1300S_20(1.55)
B = Fiber_DCF1300S_33(1.55)
C = Fiber_2028M12(1.55)

Clad = Fused2( Radius = 62.5, Fusion = 0.9, Index = Fused_silica(1.55))



Clad0 = Circle( Position = Clad.C[0], Radi = A.rClad*2, Index = A.nClad )
Core0 = Circle( Position = Clad.C[0], Radi = A.rCore, Index = A.nCore )


Clad1 = Circle( Position = Clad.C[1], Radi = B.rClad, Index = B.nClad )
Core1 = Circle( Position = Clad.C[1], Radi = B.rCore, Index = B.nCore )


SMF28 = Geometry(Objects = [Clad, Clad0, Clad1, Core0, Core1],
                 Xbound  = [-90, 90],
                 Ybound  = [-90, 90],
                 Nx      = 30,
                 Ny      = 30,
                 GConv   = 0)

#SMF28.Plot()

Sol = SuPySolver(Coupler    = SMF28,
                 Tolerance  = 1e-30,
                 MaxIter    = 10000,
                 nMode      = 5,
                 sMode      = 3,
                 Error      = 2)

SuperSet = Sol.GetModes(wavelength    = 1.55,
                          Nstep         = 300,
                          ITRi          = 1,
                          ITRf          = 0.1,
                          RightSymmetry = 0,
                          TopSymmetry   = 0,
                          Sorting       = 'Field')


class CouplerProfile(object):
    def __init__(self, L0, N, ITRf, Symmetric=True):
        self.ITRf    = ITRf
        self.L0      = L0
        self.N       = N
        self.Zend    = -L0 * np.log(ITRf)

        self.Z       = np.linspace(0,self.Zend,N)
        self.Rho     = np.exp(-(self.Z/self.L0))
        foo          = lambda x: np.exp(-(x/self.L0))

        if Symmetric:
            self.Rho     = np.concatenate((self.Rho, self.Rho[::-1]))
            self.Z       = np.concatenate((self.Z, self.Z+self.Z[-1]))

    def Plot(self):
        fig, ax = plt.subplots(1,1)
        ax.plot(Profile.Z, Profile.Rho, 'k', label='Profile')
        ax.grid()
        ax.legend()
        plt.show()

Profile = CouplerProfile(L0=1e3, N = 300, ITRf=0.1, Symmetric=False)
#Profile.Plot()

Mode = SuperSet.GetMode(Amplitude=[1/np.sqrt(2),1/np.sqrt(2),0])

Mode.Propagate(Profile=Profile.Rho, Distance=Profile.Z, atol = 1e-10, rtol=1e-6)



fig, ax = plt.subplots(1,1)
for i in range(Mode.Amplitudes.shape[0]):
    ax.plot(Mode.Z, np.abs(Mode.Amplitudes[i]), label=f"Mode {i}")

ax0 = ax.twinx()
ax0.plot(Profile.Z, Profile.Rho, 'k', label='Profile')
ax.grid()
ax0.legend()
ax.legend()
plt.show()
