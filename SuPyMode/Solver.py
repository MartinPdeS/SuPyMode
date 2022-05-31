import numpy as np

from SuPyMode.SuperMode             import SuperSet, SetSlice, SuperMode
from SuPyMode.bin.EigenSolver       import EigenSolving
from SuPyMode.Tools.utils           import Axes



class SuPySolver(object):
    """ This object corresponds to the solver.
    It solves the eigenvalues problems for a given geometry.

    """
    _Sorting = None

    def __init__(self, Coupler, Tolerance, MaxIter, Error=2,  debug='INFO'):

        Coupler.CreateMesh()
        self.Geometry     = Coupler
        self.Tolerance    = Tolerance
        self.MaxIter      = MaxIter
        self.Error        = Error
        self.CppSolvers   = []



    def InitBinding(self, Symmetries: dict, Wavelength: float, nMode: int, sMode: int):

        CppSolver = EigenSolving(Mesh       = self.Geometry.mesh,
                                 Gradient   = self.Geometry.Gradient().ravel(),
                                 nMode      = nMode,
                                 sMode      = sMode,
                                 MaxIter    = self.MaxIter,
                                 Tolerance  = self.Tolerance,
                                 dx         = self.Axes.dx,
                                 dy         = self.Axes.dy,
                                 Wavelength = Wavelength,
                                 Debug      = False
                                 )

        CppSolver.TopSymmetry    = Symmetries['Top']
        CppSolver.BottomSymmetry = Symmetries['Bottom']
        CppSolver.LeftSymmetry   = Symmetries['Left']
        CppSolver.RightSymmetry  = Symmetries['Right']


        CppSolver.ComputeLaplacian(self.Error)

        self.CppSolvers.append(CppSolver)

        return CppSolver


    def GetSuperSet(self,
                    Wavelength:     float,
                    Nstep:          int,
                    nMode:          int,
                    sMode:          int,
                    ITRi:           float = 1.0,
                    ITRf:           float = 0.1,
                    Symmetries:     dict = 0,
                    Sorting:        str   = 'Index',
                    ):

        CppSolver  = self.InitBinding(Symmetries, Wavelength, nMode, sMode)

        self.ITRList = np.linspace(ITRi, ITRf, Nstep)

        CppSolver.LoopOverITR(ITR=self.ITRList, ExtrapolationOrder=3)

        CppSolver.SortModes(Sorting=Sorting)

        Set = SuperSet(ParentSolver=self, CppSolver=CppSolver)

        self.PopulateSuperSet(Set, CppSolver)

        self.PopulateSuperModes(Set, CppSolver)

        return Set



    def PopulateSuperSet(self, Set, CppSolver):
        for m in range(CppSolver.sMode):
            Set.SuperModes.append( SuperMode(ParentSet=Set, ModeNumber=m, CppSolver=CppSolver)  )


    def PopulateSuperModes(self, Set, CppSolver):
        for SliceNumber, _ in enumerate(self.ITRList):
            Fields, Betas = CppSolver.GetSlice(SliceNumber)

            for mode, (field, beta) in enumerate( zip(Fields, Betas) ):
                Set[mode].Append(Field       = field,
                                 Index       = beta / self.Axes.k,
                                 Beta        = beta,
                                 SliceNumber = SliceNumber)


    def GetCoupling_(self):
        for Solver in self.CppSolvers:
            print(Solver)
        #Coupling = self.CppSolver.ComputingCoupling()


    @property
    def Axes(self):
        return self.Geometry.Axes



# ---
