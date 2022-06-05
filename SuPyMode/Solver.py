import numpy as np

from SuPyMode.SuperMode             import SuperSet, SuperMode
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
        #CppSolver1  = self.InitBinding({'Right': +1, 'Left': 0, 'Top': -1, 'Bottom': 0}, Wavelength, nMode, sMode)

        self.ITRList = np.linspace(ITRi, ITRf, Nstep)



        CppSolver.LoopOverITR(ITR=self.ITRList, ExtrapolationOrder=3)
        #CppSolver1.LoopOverITR(ITR=self.ITRList, ExtrapolationOrder=3)



        CppSolver.SortModes(Sorting=Sorting)
        #CppSolver1.SortModes(Sorting=Sorting)

        Set = SuperSet(ParentSolver=self)

        self.PopulateSuperSet(Set, [CppSolver])

        return Set



    def PopulateSuperSet(self, Set, CppSolvers):
        for CppSolver in CppSolvers:
            for BindingNumber in range(CppSolver.sMode):
                Set.AppendSuperMode(CppSolver=CppSolver, BindingNumber=BindingNumber)







    @property
    def Axes(self):
        return self.Geometry.Axes



# ---
