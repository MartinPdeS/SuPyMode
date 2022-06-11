import numpy as np

from SuPyMode.SuperMode             import SuperSet, SuperMode
from SuPyMode.bin.EigenSolver       import EigenSolving
from SuPyMode.Tools.utils           import Axes



class SuPySolver(object):
    """ This object corresponds to the solver.
    It solves the eigenvalues problems for a given geometry.

    """
    _Sorting = None

    def __init__(self, Coupler, Tolerance, MaxIter, Error=2,  Debug=True):

        Coupler.CreateMesh()
        self.Geometry     = Coupler
        self.Tolerance    = Tolerance
        self.MaxIter      = MaxIter
        self.Error        = Error
        self.SolverNumber = 0
        self.Debug        = Debug




    def InitBinding(self, Symmetries: dict, Wavelength: float, nMode: int, sMode: int):

        CppSolver = EigenSolving(Mesh       = self.Geometry._Mesh,
                                 Gradient   = self.Geometry.Gradient().ravel(),
                                 nMode      = nMode,
                                 sMode      = sMode,
                                 MaxIter    = self.MaxIter,
                                 Tolerance  = self.Tolerance,
                                 dx         = self.Axes.dx,
                                 dy         = self.Axes.dy,
                                 Wavelength = Wavelength,
                                 Debug      = self.Debug
                                 )

        CppSolver.TopSymmetry    = Symmetries['Top']
        CppSolver.BottomSymmetry = Symmetries['Bottom']
        CppSolver.LeftSymmetry   = Symmetries['Left']
        CppSolver.RightSymmetry  = Symmetries['Right']


        CppSolver.ComputeLaplacian(self.Error)


        return CppSolver


    def CreateSuperSet(self, Wavelength: float, NStep: int, ITRi: float, ITRf: float):
        self.Wavelength = Wavelength
        self.NStep      = NStep
        self.ITRi       = ITRi
        self.ITRf       = ITRf
        self.ITRList    = np.linspace(ITRi, ITRf, NStep)
        self.Set        = SuperSet(ParentSolver=self)


    def AddModes(self,
                    nMode:          int,
                    sMode:          int,
                    Symmetries:     dict = 0,
                    Sorting:        str   = 'Index',
                    ):

        CppSolver  = self.InitBinding(Symmetries, self.Wavelength, nMode, sMode)

        CppSolver.LoopOverITR(ITR=self.ITRList, ExtrapolationOrder=3)

        CppSolver.SortModes(Sorting='Field')

        #CppSolver.ComputeCouplingAdiabatic()


        for BindingNumber in range(CppSolver.sMode):
            self.Set.AppendSuperMode(CppSolver=CppSolver, BindingNumber=BindingNumber, SolverNumber=self.SolverNumber)

        self.SolverNumber += 1




    def GetSet(self):

        return self.Set




    @property
    def Axes(self):
        return self.Geometry.Axes



# ---
