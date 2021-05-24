import numpy               as np
import copy                as cp
import matplotlib.pyplot   as plt
from itertools             import combinations, combinations_with_replacement as Combinations
from mayavi                import mlab

from SuPyModes.Config      import *
from SuPyModes.utils       import RecomposeSymmetries, SortSuperSet
from SuPyModes.BaseClass   import SetPlots, SetProperties
from SuPyModes.Special     import Overlap, GeoGradient, ModeCoupling, ModeAdiabatic







class SuperMode(object):

    def __init__(self, Name, Profile):
        self.Name       = Name
        self.ITR        = []
        self.Index      = []
        self.Field      = []
        self.Beta       = []
        self.xSym       = []
        self.ySym       = []
        self.GeoID      = None
        self.Axes       = []
        self.Profile    = Profile
        self._Adiabatic = None
        self._Coupling  = None


    def Append(self, **kwargs):
        for property in PROPERTIES:
            getattr(self, property).append( kwargs[property] )


    @property
    def DegenerateFactor(self):
        Factor = 1

        if self.xSym[0] in [1,-1]: Factor *= 2

        if self.ySym[0] in [1,-1]: Factor *= 2

        return Factor



    def GetCoupling(self, SuperMode):
        C = []
        for n, itr in enumerate(self.ITR[:-1]):

            C.append( ModeCoupling(SuperMode0 = self,
                                   SuperMode1 = SuperMode,
                                   k          = self.Axes[n].Direct.k,
                                   Profile    = self.Profile,
                                   iter       = n) )

        return C


    def GetAdiabatic(self, SuperMode):
        A = []
        for n, itr in enumerate(self.ITR[:-1]):
            A.append( ModeAdiabatic(SuperMode0 = self,
                                    SuperMode1 = SuperMode,
                                    k          = self.Axes[n].Direct.k,
                                    Profile    = self.Profile,
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


    def FullField(self, iter):
        Field      = self.Field[iter]
        Symmetries = [self.xSym[iter], self.ySym[iter]]

        Field, xAxis, yAxis = RecomposeSymmetries(Input      = Field,
                                                  Symmetries = Symmetries,
                                                  Xaxis      = self.Axes[iter].Direct.X,
                                                  Yaxis      = self.Axes[iter].Direct.Y)

        return Field, xAxis, yAxis



class SuperSet(SetProperties, SetPlots):
    def __init__(self, IndexProfile, NSolutions, ITR):
        self.IndexProfile = IndexProfile
        self.NSolutions   = NSolutions
        self.SuperModes   = []
        self.ITR          = ITR
        self.Symmetries   = None
        self._Coupling     = None
        self._Adiabatic    = None
        self.Init()

        self.combinations = tuple(combinations( np.arange(NSolutions), 2 ) )
        self.Combinations = tuple(Combinations( np.arange(NSolutions), 2 ) )


    def Init(self):
        for solution in range(self.NSolutions):
            supermode = SuperMode(Profile = self.IndexProfile,
                                  Name = f"Mode {solution}")

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
        pass
        #self.SuperModes = SortFields(self.SuperModes)
        #for iter, _ in enumerate( self.ITR ):
        #    self.OrderingModes(iter)

    def Debug(self):
        for n, itr in enumerate( self.ITR ):
            self.Plot('Fields', iter=n)


    def OrderingModes(self, iter):
        lst = [ self[mode].Index[iter] for mode in range( self.NSolutions ) ]

        sortedIndex = sorted(range(len(lst)), key = lambda k: lst[k])

        self.SortedIndex = sortedIndex[::-1]

        temporary = cp.deepcopy( self.SuperModes )

        for n, m  in enumerate( self.SortedIndex ):
            if n != m:
                self.SwapProperties(self[n], temporary[m], iter)
            # plt.figure()
            # plt.plot(self[n].Index)
            # plt.plot(self.SuperModes[n].Index,'*')
            # plt.show()
























# -
