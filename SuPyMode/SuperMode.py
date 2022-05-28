import logging
import os
import numpy               as np
import copy                as cp
from scipy.integrate       import solve_ivp
from scipy.interpolate     import interp1d
from mayavi                import mlab

from SuPyMode.Tools.BaseClass   import SetProperties, SetPlottings
from SuPyMode.Tools.utils       import RecomposeSymmetries

Mlogger = logging.getLogger(__name__)

class SuperMode(object):

    def __init__(self, Name, ParentSet):
        self.Name       = Name
        self.Slice      = []
        self._Adiabatic = None
        self._Coupling  = None
        self.Parent     = ParentSet

    def IterateSlice(self):
        for n, slice in enumerate(self.Slice):
            yield n, slice

    @property
    def Geometry(self):
        return self.Parent.Geometry

    def Append(self, Field, Index, Beta):
        self.Slice.append(SetSlice(Field, Index=Index, Beta=Beta, ParentSet=self.Parent))


    def __getitem__(self, N):
        return self.Slice[N]


    def __setitem__(self, N, val):
        self.Slice[N] = val


    def PlotPropagation(self):
        image, _, _ = self.FullField(0)

        image = np.abs(image)

        mlab.surf(image, warp_scale="auto")

        @mlab.animate(delay=100)
        def anim_loc():
            for i in range(len(self.Field)):
                image, _, _ = self.FullField(i)
                mlab.surf( np.abs(image), warp_scale="auto")
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


    def FullField(self, iter):
        Field      = self.Slice[iter]

        return RecomposeSymmetries(Input      = Field,
                                   Symmetries = ParentSet.Symmetries,
                                   Axes       = ParentSet.Geometry.Axes)



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



class SuperSet(SetProperties, SetPlottings):
    def __init__(self, Parent):
        self.Parent        = Parent
        self.Init()


    def IterateSuperMode(self):        
        for n, supermode in enumerate(self.SuperModes):
            yield n, supermode


    def ComputeCoupling(self):
        return self.CppSolver.ComputingCoupling()


    def Init(self):
        for solution in range(self.sMode):
            supermode = SuperMode(Name=f"Mode {solution}", ParentSet=self)
            supermode.number = solution
            self.SuperModes.append(supermode)


    def ComputeM(self, CouplingFactor):
        shape = self.Beta.shape
        M     = np.zeros( [shape[0], shape[1], shape[1]] )
        for iter in range(shape[0]):
            beta = self.Beta[iter]
            M[iter] = CouplingFactor[iter] * self.Coupling[iter] + beta * np.identity(shape[1])

        return M


    def GetMode(self, Amplitude, Name = 1):
        return _SuperMode(Amplitude=Amplitude, Name=Name, ParentSet=self)


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

    @property
    def Geometry(self):
        return self.Parent.Geometry

    @property
    def Axes(self):
        return self.Parent.Geometry.Axes

    @property
    def sMode(self):
        return self.Parent.sMode

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





class SetSlice(np.ndarray):

    def __new__(cls, Field, Index, Beta, ParentSet):
        self      = Field.view(SetSlice)

        return self


    def __init__(self, Field, Index, Beta, ParentSet):

        self.Field  = Field
        self.Index  = Index
        self.Beta   = Beta
        self.ParentSet = ParentSet


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
    def Axes(self):
        return self.ParentSet.Geometry.Axes

    @property
    def Symmetries(self):
        return self.ParentSet.Symmetries

    @property
    def Norm(self):
        return (self.Field**2).sum()

    def Normalize(self):
        self.Field *= 1 / self.Norm
        self.Index *= 1 / self.Norm

    def GetField(self, Symmetries=True):
        if Symmetries:
            return RecomposeSymmetries(self.Field, self.Symmetries, self.Axes)
        else:
            self.Field, self.Axes


    def __plot__(self, ax, title=None):

        Field, xaxis, yaxis = RecomposeSymmetries(self.Field, self.Symmetries, self.Axes)


        ax.pcolormesh(yaxis, xaxis, Field.T, shading='auto', cmap='bwr')

        ax.set_ylabel(r'Y-distance [$\mu$m]', fontsize=6)

        ax.set_xlabel(r'X-distance [$\mu$m]', fontsize=6)

        ax.set_aspect('equal')
        if title:
            ax.set_title(title, fontsize=8)


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
