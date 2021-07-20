import logging
import numpy as np

from SuPyModes.SuperMode             import SuperSet, ModeSlice
from SuPyModes.includes.EigenSolver  import EigenSolving
from SuPyModes.utils                 import Axes
np.set_printoptions(precision=1, linewidth=350)
logging.basicConfig(level=logging.INFO)

Mlogger = logging.getLogger(__name__)


class SuPySolver(object):
    """ This object corresponds to the solutioner.
    It solves the eigenvalues problems for a given geometry.

    """
    def __init__(self, Coupler,  debug='INFO'):
        Mlogger.setLevel(getattr(logging, debug))
        self.Geometry     = Coupler



    def GetModes(self,
                  wavelength: float,
                  Nstep:      int   = 2,
                  Nsol:       int   = 5,
                  ITRi:       float = 1.0,
                  ITRf:       float = 0.1,
                  Ysym:       int   = 0,
                  Xsym:       int   = 0,
                  tolerance:  float = 1e-30,
                  error:      int   = 2,
                  naming:     bool  = False,
                  debug:      bool  = False):


        self.CppSolver = EigenSolving(Mesh      = self.Geometry.mesh,
                                      Gradient  = self.Geometry.Gradient().T.ravel(),
                                      nMode     = Nsol+3,
                                      sMode     = Nsol,
                                      MaxIter   = 10000,
                                      Tolerance = tolerance)

        self.CppSolver.Lambda = wavelength

        self.Geometry.Axes.wavelength = wavelength

        self.CppSolver.dx     = self.Geometry.Axes.dx

        self.CppSolver.dy     = self.Geometry.Axes.dy

        self.Geometry.Axes.Symmetries = [Xsym, Ysym]

        self.Geometry.ITRList = np.linspace(ITRi, ITRf, Nstep)

        self.Set = SuperSet(NSolutions = Nsol, Geometry = self.Geometry)

        self.Solve(self.Geometry.ITRList, Nsol)

        self.Set.CppSolver = self.CppSolver

        return self.Set



    def Solve(self, iteration_list: list, Nsol=1):

        self.CppSolver.ComputeLaplacian()

        self.CppSolver.LoopOverITR(ITR = iteration_list, ExtrapolationOrder = 3)

        self.CppSolver.SortModesFields()

        for n, _ in enumerate(iteration_list):
            Fields, Betas = self.CppSolver.GetSlice(n)

            self.AppendSlice(Betas[:Nsol], Fields[:Nsol,...].swapaxes(1,2))




    def AppendSlice(self, betas, vectors):

        for solution in range(len(betas)):

            index = betas[solution] / self.Geometry.Axes.k

            Field = vectors[solution,...]

            self.Set[solution].Slice.append(ModeSlice( Field, self.Geometry.Axes, Index=index, Beta=betas[solution] ),)


    def GetCoupling(self):
        Coupling = self.CppSolver.ComputingCoupling()


# ---
