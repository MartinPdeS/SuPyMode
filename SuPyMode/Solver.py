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



    def InitBinding(self, Symmetries: dict, Sorting: str, Wavelength: float):

        CppSolver = EigenSolving(Mesh      = self.Geometry.mesh,
                                 Gradient  = self.Geometry.Gradient().ravel(),
                                 nMode     = self.nMode,
                                 sMode     = self.sMode,
                                 MaxIter   = self.MaxIter,
                                 Tolerance = self.Tolerance,
                                 )

        CppSolver.dx             = self.Axes.dx
        CppSolver.dy             = self.Axes.dy
        CppSolver.Lambda         = Wavelength
        CppSolver.TopSymmetry    = Symmetries['Top']
        CppSolver.BottomSymmetry = Symmetries['Bottom']
        CppSolver.LeftSymmetry   = Symmetries['Left']
        CppSolver.RightSymmetry  = Symmetries['Right']

        CppSolver.ComputeLaplacian(self.Error)

        return CppSolver


    def GetModes(self,
                 Wavelength:     float,
                 Nstep:          int   = 2,
                 ITRi:           float = 1.0,
                 ITRf:           float = 0.1,
                 Symmetries:     dict = 0,
                 Sorting:        str   = 'Field'):

        CppSolver = self.InitBinding(Symmetries, Sorting, Wavelength)

        self.ITRList = np.linspace(ITRi, ITRf, Nstep)

        CppSolver.LoopOverITR(ITR=self.ITRList, ExtrapolationOrder = 3)

        self.SortModes(CppSolver=CppSolver, Sorting=Sorting)

        return self.MakeSuperSet(CppSolver)


    def SortModes(self, CppSolver, Sorting):
        if Sorting == 'Field':
            CppSolver.SortModesFields()

        elif Sorting == 'Index':
            CppSolver.SortModesIndex()


    def MakeSuperSet(self, CppSolver):

        Set = SuperSet(ParentSolver=self, CppSolver=CppSolver)

        for n, _ in enumerate(self.ITRList):

            Fields, Betas = CppSolver.GetSlice(n)

            Fields.swapaxes(1,2)

            for solution in range(self.sMode):

                index = Betas[solution] / self.Axes.k

                Field = Fields[solution,...]

                Set[solution].Append(Field       = Field,
                                     Index       = index,
                                     Beta        = Betas[solution],
                                     SliceNumber = n)

        return Set


    def GetCoupling(self):
        Coupling = self.CppSolver.ComputingCoupling()

    @property
    def Axes(self):
        return self.Geometry.Axes


    @property
    def Sorting(self):
        return self._Sorting

    @Sorting.setter
    def Sorting(self, value):
        assert value in ['Field', 'Index'], "Sorting can only be done taking account of 'Field' or 'Index'"
        self._Sorting = value


# ---
