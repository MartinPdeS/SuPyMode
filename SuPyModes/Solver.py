import logging
import numpy as np

from SuPyModes.SuperMode             import SuperSet, ModeSlice
from SuPyModes.includes.EigenSolver  import EigenSolving
from SuPyModes.utils                 import Axes

logging.basicConfig(level=logging.INFO)

Mlogger = logging.getLogger(__name__)


class SuPySolver(object):
    """ This object corresponds to the solutioner.
    It solves the eigenvalues problems for a given geometry.

    """
    def __init__(self, Coupler, Tolerance, MaxIter, nMode, sMode,  debug='INFO'):
        Mlogger.setLevel(getattr(logging, debug))
        self.Geometry     = Coupler
        self.Geometry.CreateMesh()
        self.Tolerance    = Tolerance
        self.MaxIter      = MaxIter
        self.nMode        = nMode
        self.sMode        = sMode

        self.CppSolver = EigenSolving(Mesh      = self.Geometry.mesh,
                                      Gradient  = self.Geometry.Gradient().T.ravel(),
                                      nMode     = self.nMode,
                                      sMode     = self.sMode,
                                      MaxIter   = self.MaxIter,
                                      Tolerance = self.Tolerance)

        self.CppSolver.dx     = self.Geometry.Axes.dx
        self.CppSolver.dy     = self.Geometry.Axes.dy
        self.CppSolver.ComputeLaplacian()


    def GetModes(self,
                 wavelength:     float,
                 Nstep:          int   = 2,
                 ITRi:           float = 1.0,
                 ITRf:           float = 0.1,
                 LeftSymmetry:   int   = 0,
                 RightSymmetry:  int   = 0,
                 TopSymmetry:    int   = 0,
                 BottomSymmetry: int   = 0,
                 Sorting:        str   = 'Fields'):

        assert LeftSymmetry in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"
        assert RightSymmetry in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"
        assert TopSymmetry in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"
        assert BottomSymmetry in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"

        assert Sorting in ['Fields', 'Index'], "Sorting can only be done taking account of 'Fields' or 'Index'"

        self.CppSolver.LeftSymmetry   = self.Geometry.Axes.LeftSymmetry    = LeftSymmetry
        self.CppSolver.RightSymmetry  = self.Geometry.Axes.RightSymmetry   = RightSymmetry
        self.CppSolver.TopSymmetry    = self.Geometry.Axes.TopSymmetry     = TopSymmetry
        self.CppSolver.BottomSymmetry = self.Geometry.Axes.BottomSymmetry  = BottomSymmetry

        self.CppSolver.Lambda = wavelength

        self.Geometry.Axes.wavelength = wavelength

        self.Geometry.Axes.Symmetries = [0, 0]

        self.Geometry.ITRList = np.linspace(ITRi, ITRf, Nstep)

        self.Set = SuperSet(NSolutions = self.sMode, Geometry = self.Geometry)

        self.Solve(self.Geometry.ITRList, self.sMode)

        self.Set.CppSolver = self.CppSolver

        if Sorting == 'Fields':
            self.CppSolver.SortModesFields()

        elif Sorting == 'Index':
            self.CppSolver.SortModesIndex()

        return self.Set


    def Solve(self, iteration_list: list, Nsol=1):

        self.CppSolver.LoopOverITR(ITR = iteration_list, ExtrapolationOrder = 1)

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





    @property
    def TopSymmetry(self):
        return self.CppSolver.TopSymmetry

    @TopSymmetry.setter
    def TopSymmetry(self, value):
        self.CppSolver.TopSymmetry = value


    @property
    def BottomSymmetry(self):
        return self.CppSolver.BottomSymmetry

    @BottomSymmetry.setter
    def BottomSymmetry(self, value):
        self.CppSolver.BottomSymmetry = value


    @property
    def LeftSymmetry(self):
        return self.CppSolver.LeftSymmetry

    @LeftSymmetry.setter
    def LeftSymmetry(self, value):
        self.CppSolver.LeftSymmetry = value


    @property
    def RightSymmetry(self):
        return self.CppSolver.RightSymmetry

    @RightSymmetry.setter
    def RightSymmetry(self, value):
        self.CppSolver.RightSymmetry = value

# ---
