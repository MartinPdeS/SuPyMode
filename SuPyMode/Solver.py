import numpy as np

from SuPyMode.SuperMode             import SuperSet, SetSlice
from SuPyMode.bin.EigenSolver       import EigenSolving
from SuPyMode.Tools.utils           import Axes



class SuPySolver(object):
    """ This object corresponds to the solver.
    It solves the eigenvalues problems for a given geometry.

    """
    _Sorting = None

    def __init__(self, Coupler, Tolerance, MaxIter, nMode, sMode, Error=2,  debug='INFO', Debug=False):

        Coupler.CreateMesh()
        self.Geometry     = Coupler
        self.Tolerance    = Tolerance
        self.MaxIter      = MaxIter
        self.nMode        = nMode
        self.sMode        = sMode
        self.Error        = Error
        self.Debug        = Debug
        self.InitBinding()


    def InitBinding(self):
        self.CppSolver = EigenSolving(Mesh      = self.Geometry.mesh,
                                      Gradient  = self.Geometry.Gradient().ravel(),
                                      nMode     = self.nMode,
                                      sMode     = self.sMode,
                                      MaxIter   = self.MaxIter,
                                      Tolerance = self.Tolerance,
                                      Debug     = self.Debug)

        self.CppSolver.dx     = self.Axes.dx
        self.CppSolver.dy     = self.Axes.dy
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

        self.LeftSymmetry   = LeftSymmetry
        self.RightSymmetry  = RightSymmetry
        self.TopSymmetry    = TopSymmetry
        self.BottomSymmetry = BottomSymmetry

        self.Sorting        = Sorting

        self.Geometry.Symmetries = self.Symmetries

        self.Wavelength = wavelength

        self.ITRList = np.linspace(ITRi, ITRf, Nstep)

        self.CppSolver.LoopOverITR(ITR = self.Geometry.ITRList, ExtrapolationOrder = 3)

        return self.MakeSuperSet()



    def SortModes(self, Sorting):
        if Sorting == 'Field':
            self.CppSolver.SortModesFields()

        elif Sorting == 'Index':
            self.CppSolver.SortModesIndex()


    def MakeSuperSet(self):

        Set = SuperSet(ParentSolver=self)

        self.SortModes(self.Sorting)

        for n, _ in enumerate(self.Geometry.ITRList):

            Fields, Betas = self.CppSolver.GetSlice(n)

            Fields.swapaxes(1,2)

            for solution in range(self.sMode):

                index = Betas[solution] / self.Axes.k

                Field = Fields[solution,...]

                Set[solution].Append(Field       = Field,
                                     Index       = index,
                                     Beta        = Betas[solution],
                                     SliceNumber = n)

                Set[solution].AddSymmetries(Left=self.LeftSymmetry, Right=self.RightSymmetry, Top=self.TopSymmetry, Bottom=self.BottomSymmetry)

        return Set


    def GetCoupling(self):
        Coupling = self.CppSolver.ComputingCoupling()

    @property
    def Axes(self):
        return self.Geometry.Axes

    @property
    def ITRList(self):
        return self.Geometry.ITRList

    @ITRList.setter
    def ITRList(self, value):
        self.Geometry.ITRList = value

    @property
    def Wavelength(self):
        return self.CppSolver.Lambda

    @Wavelength.setter
    def Wavelength(self, value):
        self.Geometry.Axes.Wavelength = value
        self.CppSolver.Lambda = value

    @property
    def Sorting(self):
        return self._Sorting

    @Sorting.setter
    def Sorting(self, value):
        assert value in ['Field', 'Index'], "Sorting can only be done taking account of 'Field' or 'Index'"
        self._Sorting = value

    @property
    def Symmetries(self):
        return {'Left'   : self.LeftSymmetry,
                'Right'  : self.RightSymmetry,
                'Top'    : self.TopSymmetry,
                'Bottom' : self.BottomSymmetry}

    @property
    def TopSymmetry(self):
        return self.CppSolver.TopSymmetry

    @TopSymmetry.setter
    def TopSymmetry(self, value):
        assert value in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"
        self.CppSolver.TopSymmetry = value

    @property
    def BottomSymmetry(self):
        return self.CppSolver.BottomSymmetry

    @BottomSymmetry.setter
    def BottomSymmetry(self, value):
        assert value in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"
        self.CppSolver.BottomSymmetry = value

    @property
    def LeftSymmetry(self):
        return self.CppSolver.LeftSymmetry

    @LeftSymmetry.setter
    def LeftSymmetry(self, value):
        assert value in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"
        self.CppSolver.LeftSymmetry = value

    @property
    def RightSymmetry(self):
        return self.CppSolver.RightSymmetry

    @RightSymmetry.setter
    def RightSymmetry(self, value):
        assert value in [-1,0,1], "Symmetries can only take the following values -1 [antisymmetric], 0 [no symmetries], 1 [symmetric]"
        self.CppSolver.RightSymmetry = value

# ---
