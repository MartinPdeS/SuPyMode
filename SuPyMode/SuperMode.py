import logging
import os
import numpy               as np
import copy                as cp
from scipy.integrate       import solve_ivp
from scipy.interpolate     import interp1d
from mayavi                import mlab

from SuPyMode.Tools.BaseClass   import SetProperties, SetPlottings


Mlogger = logging.getLogger(__name__)


class SuperSet(SetProperties, SetPlottings):
    def __init__(self, ParentSolver):
        self.ParentSolver   = ParentSolver
        self.SuperModes     = [SuperMode(ParentSet=self, ModeNumber=m) for m in range(self.ParentSolver.sMode)]


    def IterateSuperMode(self):
        for n, supermode in enumerate(self.SuperModes):
            yield supermode


    def ComputeCoupling(self):
        return self.CppSolver.ComputingCoupling()


    def ComputeM(self, CouplingFactor):
        shape = self.Beta.shape
        M     = np.zeros( [shape[0], shape[1], shape[1]] )
        for iter in range(shape[0]):
            beta = self.Beta[iter]
            M[iter] = CouplingFactor[iter] * self.Coupling[iter] + beta * np.identity(shape[1])

        return M


    def ComputeCouplingFactor(self, Length):
        dx =  Length/(self.Geometry.ITRList.size)

        dITR = np.gradient(np.log(self.Geometry.ITRList), 1)

        factor = dITR/dx

        return factor/10


    def Propagate(self, Amplitude, Length, **kwargs):
        Amplitude = np.asarray(Amplitude)

        Distance = np.linspace(0, Length, self.Geometry.ITRList.size)

        Factor = self.ComputeCouplingFactor(Length)


        M = self.ComputeM(CouplingFactor=Factor)

        Minterp = interp1d(Distance, M, axis=0)

        def foo(t, y):
            return 1j * Minterp(t).dot(y)

        sol = solve_ivp(foo,
                        y0       = Amplitude.astype(complex),
                        t_span   = [0, Length],
                        method   = 'RK45',
                        **kwargs)

        return sol


    def PlotPropagation(self, Modes):
        Mode0 = self[0]
        Mode1 = self[1]

        Field = Mode0[0].GetFullField(Axes=False) + Mode1[0].GetFullField(Axes=False)

        surface = mlab.surf( np.abs( Field ) , warp_scale="auto" )

        @mlab.animate(delay=100)
        def anim_loc():
            for n, _ in Mode0.IterateSlice():
                Field = Mode0[n].GetFullField(Axes=False) + Mode1[n].GetFullField(Axes=False)

                surface.mlab_source.scalars = np.abs( Field )

                yield

        anim_loc()
        mlab.show
    

    @property
    def Size(self):
        return len(self.SuperModes)

    @property
    def Geometry(self):
        return self.ParentSolver.Geometry

    @property
    def ITRList(self):
        return self.ParentSolver.Geometry.ITRList

    @property
    def Axes(self):
        return self.ParentSolver.Geometry.Axes


    def GetCoupling(self, Mode0: int, Mode1: int):
        return self.Coupling[:, Mode0, Mode1]


    def GetAdiabatic(self, Mode0: int, Mode1: int):
        return self.Adiabatic[:, Mode0, Mode1]


    def Append(self, Object):
        self.SuperModes.append(Object)


    def __getitem__(self, N):
        return self.SuperModes[N]


    def __setitem__(self, N, val):
        self.SuperModes[N] = val


    def __copy__(self):
        to_be_copied = ['SuperModes']

        copy_ = SuperSet(Parent)

        for attr in self.__dict__:

            if attr in to_be_copied:
                copy_.__dict__[attr] =  cp.deepcopy( self.__dict__[attr] )
            else:
                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_


















class SuperMode(object):

    def __init__(self, ParentSet, ModeNumber):
        self.Name       = f"Mode: {ModeNumber}"
        self.ModeNumber = ModeNumber
        self.Slice      = []
        self.ParentSet  = ParentSet
        self.Index      = self.ParentSet.Index[:, self.ModeNumber]
        self.Beta       = self.ParentSet.Beta[:, self.ModeNumber]



    def AddSymmetries(self, Left, Right, Top, Bottom):
        self.LeftSymmetry   = Left
        self.RightSymmetry  = Right
        self.TopSymmetry    = Top
        self.BottomSymmetry = Bottom


    def IterateSlice(self):
        for n, slice in enumerate(self.Slice):
            yield n, slice


    @property
    def Geometry(self):
        return self.ParentSet.Geometry


    @property
    def ITRList(self):
        return self.ParentSet.Geometry.ITRList


    @property
    def Axes(self):
        return self.ParentSet.Geometry.Axes


    @property
    def Symmetries(self):
        return np.array(self.LeftSymmetry, self.RightSymmetry, self.TopSymmetry, self.BottomSymmetry)


    def CompareSymmetries(self, Other):
        assert isinstance(Other, SuperMode), "Can only compare SuperMode instance with another"
        return np.all(self.Symmetries == Other.Symmetries)


    def Append(self, Field, Index, Beta, SliceNumber):
        self.Slice.append( SetSlice(Field, Index=Index, Beta=Beta, ParentMode=self, SliceNumber=SliceNumber) )

    def GetSlices(self, SlicesNumber: list):
        return [self[slice] for slice in SlicesNumber]


    def __getitem__(self, N):
        return self.Slice[N]


    def __setitem__(self, N, val):
        self.Slice[N] = val


    def PlotPropagation(self):
        image = np.abs( self[0].GetFullField(Axes=False) )

        surface = mlab.surf(image, warp_scale="auto")

        @mlab.animate(delay=100)
        def anim_loc():
            for n, slice in self.IterateSlice():
                surface.mlab_source.scalars = np.abs( slice.GetFullField(Axes=False) )

                yield

        anim_loc()
        mlab.show()

        """
        import os
        fps = 20
        prefix = 'ani'
        ext = '.png'

        import subprocess
        animate_plots(base_directory='yolo', fname_prefix='dasda')"""



    def __copy__(self):
        to_be_copied = ['Slice']

        copy_ = SuperSet( self.Name, self.Geometry)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_


    def __deepcopy__(self, memo):
        to_be_copied = ['Slice']

        copy_ = SuperMode(Name=self.Name)

        for attr in self.__dict__:

            if attr in to_be_copied:
                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])

            else:
                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_























class SetSlice(np.ndarray):

    def __new__(cls, Field, Index, Beta, ParentMode, SliceNumber):
        self      = Field.view(SetSlice)
        return self


    def CompareSymmetries(self, Other):
        assert isinstance(Other, SetSlice), "Can only compare SetSlice instance with another"
        return np.all(self.Symmetries == Other.Symmetries)

    def __init__(self, Field, Index, Beta, ParentMode, SliceNumber):
        self.Field  = Field
        self.Index  = Index
        self.Beta   = Beta
        self.SliceNumber = SliceNumber
        self.ParentMode = ParentMode


    def __array_finalize__(self, viewed):
        pass


    def __pow__(self, other):
        assert isinstance(other, SetSlice), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    def Overlap(self, other):
        assert isinstance(other, SetSlice), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    @property
    def LeftSymmetry(self):
        return self.ParentMode.LeftSymmetry


    @property
    def ITR(self):
        return self.ParentMode.ITRList[self.SliceNumber]

    @property
    def RightSymmetry(self):
        return self.ParentMode.RightSymmetry

    @property
    def TopSymmetry(self):
        return self.ParentMode.TopSymmetry

    @property
    def BottomSymmetry(self):
        return self.ParentMode.BottomSymmetry

    @property
    def Axes(self):
        return self.ParentMode.Geometry.Axes

    @property
    def Symmetries(self):
        return self.ParentMode.Symmetries

    @property
    def Norm(self):
        return (self.Field**2).sum()

    def Normalize(self):
        self.Field *= 1 / self.Norm
        self.Index *= 1 / self.Norm

    def GetField(self, Symmetries=True):
        return self.Field, self.Axes.X, self.Axes.Y


    def GetFullField(self, Axes: bool = True):
        Xaxis = self.Axes.Y.copy()
        Yaxis = self.Axes.X.copy()
        Field = self.Field.copy()

        if self.BottomSymmetry == 1:
            Field = np.concatenate((Field[:,::-1], Field), axis=1)
            Xaxis = np.concatenate( [Xaxis + Xaxis[0] - Xaxis[-1], Xaxis] )


        elif self.BottomSymmetry == -1:
            Field = np.concatenate((-Field[:,::-1], Field), axis=1)
            Xaxis = np.concatenate( [Xaxis + Xaxis[0] - Xaxis[-1], Xaxis] )


        if self.TopSymmetry == 1:
            Field = np.concatenate((Field, Field[:,::-1]), axis=1)
            Xaxis = np.concatenate( [Xaxis, Xaxis - Xaxis[0] + Xaxis[-1]] )

        elif self.TopSymmetry == -1:
            Field = np.concatenate((Field, -Field[:,::-1]), axis=1)
            Xaxis = np.concatenate( [Xaxis, Xaxis - Xaxis[0] + Xaxis[-1]] )


        if self.RightSymmetry == 1:
            Field = np.concatenate((Field, Field[::-1,:]), axis=0)
            Yaxis = np.concatenate( [Yaxis, Yaxis - Yaxis[0] + Yaxis[-1]] )

        elif self.RightSymmetry == -1:
            Field = np.concatenate((Field, -Field[::-1,:]), axis=0)
            Yaxis = np.concatenate( [Yaxis, Yaxis - Yaxis[0] + Yaxis[-1]] )


        if self.LeftSymmetry == 1:
            Field = np.concatenate((Field[::-1,:], Field), axis=0)
            Yaxis = np.concatenate( [Yaxis + Yaxis[0] - Yaxis[-1], Yaxis] )

        elif self.LeftSymmetry == -1:
            Field = np.concatenate((-Field[::-1,:], Field), axis=0)
            Yaxis = np.concatenate( [Yaxis + Yaxis[0] - Yaxis[-1], Yaxis] )

        if Axes is True:
            return Field, Xaxis, Yaxis

        if Axes is False:
            return Field



    def __copy__(self):
        to_be_copied = ['Field', 'Index', 'Axes', 'Beta']

        copy_ = SetSlice(self, self.Axes, self.Index, self.Beta)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_



    def __deepcopy__(self, memo):
        to_be_copied = ['Field', 'Index', 'Axes', 'Beta']

        copy_ = SetSlice(self, self.Axes, self.Index, self.Beta)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_


















# -
