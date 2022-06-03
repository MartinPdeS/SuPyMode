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
        self.SuperModes     = []
        self._NextMode      = 0


    @property
    def NextMode(self):
        self._NextMode += 1
        return self._NextMode - 1

    def IterateSuperMode(self):
        for n, supermode in enumerate(self.SuperModes):
            yield supermode


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
        return self.ParentSolver.ITRList

    @property
    def Axes(self):
        return self.ParentSolver.Geometry.Axes


    def AppendSuperMode(self, CppSolver, BindingNumber):
        superMode = SuperMode(ParentSet=self, ModeNumber=self.NextMode, CppSolver=CppSolver, BindingNumber=BindingNumber )

        self.SuperModes.append( superMode )


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

    def __init__(self, ParentSet, ModeNumber, CppSolver,  BindingNumber):
        self.Name           = f"Mode: {ModeNumber}"
        self.ModeNumber     = ModeNumber
        self.BindingNumber  = BindingNumber
        self.CppSolver      = CppSolver
        self.Slice          = []
        self.ParentSet      = ParentSet


    def GetCoupling(self, OtherMode):
        return self.CppSolver.ComputingCoupling()


    def GetAdiabatic(self, OtherMode):
        if self.CppSolver != OtherMode.CppSolver:
            return np.zeros(self.CppSolver.ITRLength)

        return self.CppSolver.ComputingAdiabatic()[:, self.BindingNumber, OtherMode.BindingNumber]


    def GetIndex(self):
        return self.CppSolver.GetIndices()[:, self.BindingNumber]


    def GetBetas(self):
        return self.CppSolver.GetBetas()[:, self.BindingNumber]


    @property
    def LeftSymmetry(self):
        return self.CppSolver.LeftSymmetry

    @property
    def RightSymmetry(self):
        return self.CppSolver.RightSymmetry

    @property
    def TopSymmetry(self):
        return self.CppSolver.TopSymmetry

    @property
    def BottomSymmetry(self):
        return self.CppSolver.BottomSymmetry


    def IterateSlice(self):
        for n, slice in enumerate(self.Slice):
            yield n, slice

    @property
    def Size(self):
        return len(self.Slice)

    @property
    def Geometry(self):
        return self.ParentSet.Geometry


    @property
    def ITRList(self):
        return self.ParentSet.ParentSolver.ITRList


    @property
    def Axes(self):
        return self.ParentSet.Geometry.Axes


    @property
    def Symmetries(self):
        return np.array(self.LeftSymmetry, self.RightSymmetry, self.TopSymmetry, self.BottomSymmetry)


    def CompareSymmetries(self, Other):
        assert isinstance(Other, SuperMode), "Can only compare SuperMode instance with another"
        return np.all(self.Symmetries == Other.Symmetries)


    def AppendSlice(self, SliceNumber):
        self.Slice.append( SetSlice(ParentMode=self, SliceNumber=SliceNumber) )

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


    def GetIndex(self):
        return self.CppSolver.GetIndices()[:, self.BindingNumber]


    def GetBetas(self):
        return self.CppSolver.GetBetas()[:, self.BindingNumber]



class SetSlice():

    def CompareSymmetries(self, Other):
        assert isinstance(Other, SetSlice), "Can only compare SetSlice instance with another"
        return np.all(self.Symmetries == Other.Symmetries)


    def __init__(self, ParentMode, SliceNumber):
        self.ParentMode = ParentMode
        self.SliceNumber = SliceNumber


    @property
    def Index(self):
        return self.ParentMode.GetIndex()[self.SliceNumber]


    @property
    def Beta(self):
        return self.ParentMode.GetBeta()[self.SliceNumber]


    @property
    def Field(self):
        return self.ParentMode.CppSolver.GetFields(self.SliceNumber)[self.ParentMode.BindingNumber]


    def __pow__(self, other):
        assert isinstance(other, SetSlice), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    def Overlap(self, other):
        assert isinstance(other, SetSlice), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    @property
    def ITR(self):
        return self.ParentMode.ITRList[self.SliceNumber]


    @property
    def LeftSymmetry(self):
        return self.ParentMode.LeftSymmetry

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
