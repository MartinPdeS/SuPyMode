import numpy               as np
import copy                as cp
import matplotlib.pyplot   as plt
from itertools             import combinations, combinations_with_replacement as Combinations
from mayavi                import mlab

from SuPyModes.Config      import *
from SuPyModes.utils       import RecomposeSymmetries
from SuPyModes.BaseClass   import SetPlots, SetProperties
from SuPyModes.Special     import ModeCoupling, ModeAdiabatic



class SuperMode(object):

    def __init__(self, Name, Geometry):
        self.Name       = Name
        self.Slice      = []
        self._Adiabatic = None
        self._Coupling  = None
        self.Geometry   = Geometry



    def Append(self, **kwargs):
        self.Field.append(kwargs['Field'])


    def __getitem__(self, N):
        return self.Slice[N]


    def __setitem__(self, N, val):
        self.Slice[N] = val


    @property
    def DegenerateFactor(self):
        Factor = 1

        if self.Geometry.Axes.Symmetries[0] in [1,-1]: Factor *= 2

        if self.Geometry.Axes.Symmetries[1] in [1,-1]: Factor *= 2

        return Factor


    def GetCoupling(self, SuperMode):
        C = []
        for n, itr in enumerate(self.Geometry.ITRList[:-1]):
            C.append( ModeCoupling(SuperMode0 = self,
                                   SuperMode1 = SuperMode,
                                   k          = self.Geometry.Axes.Direct.k,
                                   Profile    = self.Geometry.mesh,
                                   iter       = n) )

        return C


    def GetAdiabatic(self, SuperMode):
        A = []
        for n, itr in enumerate(self.Geometry.ITRList[:-1]):
            A.append( ModeAdiabatic(SuperMode0 = self,
                                    SuperMode1 = SuperMode,
                                    k          = self.Geometry.Axes.Direct.k,
                                    Profile    = self.Geometry.mesh,
                                    iter       = n) )

        return A


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

        import os
        fps = 20
        prefix = 'ani'
        ext = '.png'

        import subprocess
        animate_plots(base_directory='yolo', fname_prefix='dasda')


    def FullField(self, iter):
        Field      = self.Field[iter]

        Field, xAxes, yAxes = RecomposeSymmetries(Input = Field, Axes = self.Axes)

        return Field, xAxes, yAxes


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

        copy_ = SuperMode( self.Name, self.Geometry)

        for attr in self.__dict__:

            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_



class SuperSet(SetProperties, SetPlots):
    def __init__(self, NSolutions, Geometry):
        self.NSolutions    = NSolutions
        self.SuperModes    = []
        self._Coupling     = None
        self._Adiabatic    = None
        self.Geometry      = Geometry
        self.Init()

        self.combinations = tuple(combinations( np.arange(NSolutions), 2 ) )
        self.Combinations = tuple(Combinations( np.arange(NSolutions), 2 ) )


    def Init(self):

        for solution in range(self.NSolutions):
            supermode = SuperMode(Name     = f"Mode {solution}",
                                  Geometry = self.Geometry)

            self.SuperModes.append(supermode)


    def __getitem__(self, N):
        return self.SuperModes[N]


    def __setitem__(self, N, val):
        self.SuperModes[N] = val


    def SwapProperties(self, SuperMode0, SuperMode1, N):
        S0, S1 = SuperMode0, SuperMode1

        for p in PROPERTIES:
            getattr(S0, p)[N] = getattr(S1, p)[N]


    def Ordering(self):
        for iter, _ in enumerate( self.Geometry.ITR ):
            self.OrderingModes(iter)


    def Debug(self):
        for n, itr in enumerate( self.Geometry.ITR ):
            self.Plot('Fields', iter=n)


    def OrderingModes(self, iter):
        lst = [ self[mode].Index[iter] for mode in range( self.NSolutions ) ]

        sortedIndex = sorted(range(len(lst)), key = lambda k: lst[k])

        self.SortedIndex = sortedIndex[::-1]

        temporary = cp.deepcopy( self.SuperModes )

        for n, m  in enumerate( self.SortedIndex ):
            if n != m:
                self.SwapProperties(self[n], temporary[m], iter)


    def __copy__(self):
        to_be_copied = ['SuperModes']

        copy_ = SuperSet(self.NSolutions, self.Geometry)

        for attr in self.__dict__:

            if attr in to_be_copied:
                copy_.__dict__[attr] =  cp.deepcopy( self.__dict__[attr] )
            else:
                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_






class ModeSlice(np.ndarray):

    def __new__(cls, Field, Axes, Index, Beta):
        self      = Field.view(ModeSlice)

        return self


    def __init__(self, Field, Axes, Index, Beta):
        self.Axes  = Axes
        self.Index = Index
        self.Beta  = Beta


    def __array_finalize__(self, viewed):
        pass


    def __pow__(self, other):
        assert isinstance(other, ModeSlice), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    def Overlap(self, other):
        assert isinstance(other, ModeSlice), f'Cannot multiply supermodes with {other.__class__}'

        overlap = np.abs( np.sum( np.multiply( self, other ) ) )

        return float( overlap )


    def __plot__(self, ax, title=None):
        Field, xaxis, yaxis = RecomposeSymmetries(self, self.Axes)

        ax.pcolormesh(yaxis, xaxis, Field, shading='auto')

        ax.set_ylabel(r'Y-distance [$\mu$m]', fontsize=6)

        ax.set_xlabel(r'X-distance [$\mu$m]', fontsize=6)

        ax.set_aspect('equal')
        if title:
            ax.set_title(title, fontsize=8)


    def __copy__(self):
        to_be_copied = ['Index', 'Axes', 'Beta']

        copy_ = ModeSlice(self, self.Axes, self.Index, self.Beta)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_



    def __deepcopy__(self, memo):
        to_be_copied = ['Index', 'Axes', 'Beta']

        copy_ = ModeSlice(self, self.Axes, self.Index, self.Beta)

        for attr in self.__dict__:
            if attr in to_be_copied:

                copy_.__dict__[attr] = cp.copy(self.__dict__[attr])
            else:

                copy_.__dict__[attr] = self.__dict__[attr]

        return copy_


















# -
