import numpy as np

from SuPyMode.SuperMode             import SuperSet, SetSlice
from SuPyMode.bin.EigenSolver       import EigenSolving
from SuPyMode.utils                 import Axes



class SuPySolver(object):
    """ This object corresponds to the solutioner.
    It solves the eigenvalues problems for a given geometry.

    """
    def __init__(self, Coupler, Tolerance, MaxIter, nMode, sMode, Error=2,  debug='INFO', Debug=False):

        self.Geometry     = Coupler
        self.Geometry.CreateMesh()
        self.Tolerance    = Tolerance
        self.MaxIter      = MaxIter
        self.nMode        = nMode
        self.sMode        = sMode
        self.Error        = Error
        self.Debug        = Debug

        self.CppSolver = EigenSolving(Mesh      = self.Geometry.mesh,
                                      Gradient  = self.Geometry.Gradient().ravel(),
                                      nMode     = self.nMode,
                                      sMode     = self.sMode,
                                      MaxIter   = self.MaxIter,
                                      Tolerance = self.Tolerance,
                                      Debug     = self.Debug)

        self.CppSolver.dx     = self.Geometry.Axes.dx
        self.CppSolver.dy     = self.Geometry.Axes.dy
        self.CppSolver.ComputeLaplacian(self.Error)


    def GetModes(self,
                 wavelength:     float,
                 Nstep:          int   = 2,
                 ITRi:           float = 1.0,
                 ITRf:           float = 0.1,
                 LeftSymmetry:   int   = 0,
                 RightSymmetry:  int   = 0,
                 TopSymmetry:    int   = 0,
                 BottomSymmetry: int   = 0,
                 Sorting:        str   = 'Field'):

        Symmetries = [LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry]
        for Sym in Symmetries:
            assert Sym in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"

        assert Sorting in ['Field', 'Index'], "Sorting can only be done taking account of 'Fields' or 'Index'"

        Set = SuperSet(NSolutions = self.sMode, Geometry = self.Geometry)

        self.CppSolver.LeftSymmetry   = LeftSymmetry
        self.CppSolver.RightSymmetry  = RightSymmetry
        self.CppSolver.TopSymmetry    = TopSymmetry
        self.CppSolver.BottomSymmetry = BottomSymmetry

        self.Geometry.Symmetries = {'LeftSymmetry'   : LeftSymmetry,
                                    'RightSymmetry'  : RightSymmetry,
                                    'TopSymmetry'    : TopSymmetry,
                                    'BottomSymmetry' : BottomSymmetry}

        self.CppSolver.Lambda = wavelength

        self.Geometry.Axes.wavelength = wavelength

        self.Geometry.ITRList = np.linspace(ITRi, ITRf, Nstep)

        self.CppSolver.LoopOverITR(ITR = self.Geometry.ITRList, ExtrapolationOrder = 3)

        return self.MakeSet(self.Geometry.ITRList, self.sMode, Sorting)



    def SortModes(self, Sorting):
        if Sorting == 'Field':
            self.CppSolver.SortModesFields()

        elif Sorting == 'Index':
            self.CppSolver.SortModesIndex()


    def MakeSet(self, iteration_list: list, Nsol=1, Sorting='Field'):

        Set = SuperSet(NSolutions = self.sMode, Geometry = self.Geometry)

        Set.Symmetries =  self.Geometry.Symmetries

        Set.CppSolver = self.CppSolver

        self.SortModes(Sorting)

        for n, _ in enumerate(iteration_list):

            Fields, Betas = self.CppSolver.GetSlice(n)

            Fields.swapaxes(1,2)

            for solution in range(Nsol):

                index = Betas[solution] / self.Geometry.Axes.k

                Field = Fields[solution,...]

                Set[solution].Slice.append(SetSlice( Field, Index=index, Beta=Betas[solution], ParentSet=Set ))

        return Set


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
