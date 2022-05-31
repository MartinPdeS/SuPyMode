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



    def InitBinding(self, LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry):
        CppSolver = EigenSolving(Mesh      = self.Geometry.mesh,
                                 Gradient  = self.Geometry.Gradient().ravel(),
                                 nMode     = self.nMode,
                                 sMode     = self.sMode,
                                 MaxIter   = self.MaxIter,
                                 Tolerance = self.Tolerance)

        CppSolver.dx     = self.Axes.dx
        CppSolver.dy     = self.Axes.dy
        CppSolver.ComputeLaplacian(self.Error)

        self.SetBindingSymmetries(CppSolver = CppSolver,
                                  Left      = LeftSymmetry,
                                  Right     = RightSymmetry,
                                  Top       = TopSymmetry,
                                  Bottom    = BottomSymmetry)

        return CppSolver


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

        self.CppSolver = self.InitBinding(LeftSymmetry, RightSymmetry, TopSymmetry, BottomSymmetry)
        self.CppSolver1 = self.InitBinding(+1, RightSymmetry, TopSymmetry, BottomSymmetry)

        self.Sorting        = Sorting

        self.Wavelength = wavelength

        self.ITRList = np.linspace(ITRi, ITRf, Nstep)

        self.CppSolver.LoopOverITR(ITR=self.Geometry.ITRList, ExtrapolationOrder = 3)



        return self.MakeSuperSet()


    def SetBindingSymmetries(self, CppSolver, Left, Right, Top, Bottom):
        CppSolver.TopSymmetry    = Top
        CppSolver.BottomSymmetry = Bottom
        CppSolver.LeftSymmetry   = Left
        CppSolver.RightSymmetry  = Right


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

                Set[solution].AddSymmetries(Left   = self.CppSolver.LeftSymmetry,
                                            Right  = self.CppSolver.RightSymmetry,
                                            Top    = self.CppSolver.TopSymmetry,
                                            Bottom = self.CppSolver.BottomSymmetry)

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


# ---
