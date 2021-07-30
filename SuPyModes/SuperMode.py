import logging
import numpy               as np
import copy                as cp
from itertools             import combinations, combinations_with_replacement as Combinations
from mayavi                import mlab

from SuPyModes.Config      import *
from SuPyModes.utils       import RecomposeSymmetries, Enumerate
from SuPyModes.BaseClass   import SetPlots, SetProperties



Mlogger = logging.getLogger(__name__)

class SuperMode(object):

    def __init__(self, Name, Geometry, ParentSet):
        self.Name       = Name
        self.Slice      = []
        self._Adiabatic = None
        self._Coupling  = None
        self.Geometry   = Geometry
        self.Parent     = ParentSet



    def Append(self, **kwargs):
        self.Slice.append(kwargs['Field'])


    def __getitem__(self, N):
        return self.Slice[N]


    def __setitem__(self, N, val):
        self.Slice[N] = val


    def GetCoupling(self, SuperMode):
        return self.Parent.Coupling[:, self.number, SuperMode.number]


    def GetAdiabatic(self, SuperMode):
        return self.Parent.Adiabatic[:, self.number, SuperMode.number]


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

        Field, xAxes, yAxes = RecomposeSymmetries(Input      = Field,
                                                  Symmetries = ParentSet.Symmetries,
                                                  Axes       = ParentSet.Geometry.Axes)

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
    def __init__(self, NSolutions, Geometry, debug='INFO'):
        Mlogger.setLevel(getattr(logging, debug))

        self.NSolutions    = NSolutions
        self.SuperModes    = []
        self._Coupling     = None
        self._Adiabatic    = None
        self._Index        = None
        self._Beta         = None
        self.Geometry      = Geometry
        self._M            = None

        self.Init()

        self.combinations = tuple(combinations( np.arange(NSolutions), 2 ) )
        self.Combinations = tuple(Combinations( np.arange(NSolutions), 2 ) )


    def ComputeCoupling(self):
        return self.CppSolver.ComputingCoupling()


    def Init(self):
        for solution in range(self.NSolutions):
            supermode = SuperMode(Name = f"Mode {solution}", Geometry = self.Geometry, ParentSet=self)
            supermode.number = solution
            self.SuperModes.append(supermode)


    def __getitem__(self, N):
        return self.SuperModes[N]


    def __setitem__(self, N, val):
        self.SuperModes[N] = val


    def __copy__(self):
        to_be_copied = ['SuperModes']

        copy_ = SuperSet(self.NSolutions, self.Geometry)

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
    def Norm(self):
        return (self.Field**2).sum()

    def Normalize(self):
        norm = self.Norm
        self.Field *= 1 / norm
        self.Index *= 1 / norm

    def __plot__(self, ax, title=None):

        Field, xaxis, yaxis = RecomposeSymmetries(self.Field, self.ParentSet.Symmetries, self.ParentSet.Geometry.Axes)


        ax.pcolormesh(yaxis, xaxis, Field.T, shading='auto', cmap='bwr')
        #ax.pcolormesh(Field.T, shading='auto', cmap='bwr')
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
