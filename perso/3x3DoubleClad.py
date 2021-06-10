from SuPyModes.Geometry          import Geometry, Circle, Fused3
from SuPyModes.Solver            import SuPySolver
from SuPyModes.sellmeier         import Fused_silica
from SuPyModes.utils             import NA2nCore

nCore = NA2nCore( 0.12, Fused_silica(1.55) )
nClad = NA2nCore( 0.11, Fused_silica(1.55) )
FiberA = {'Core':(nCore, 4.6), 'Clad': (nClad,19.9/2)}


nCore = NA2nCore( 0.125, Fused_silica(1.55) )
nClad = NA2nCore( 0.11, Fused_silica(1.55)  )
FiberB = {'Core':(nCore, 4.5), 'Clad': (nClad, 33/2)}



Clad = Fused3(Radius =  62.5, Fusion  = 0.8, Index   = Fused_silica(1.55))


Clad0 = Circle( Position = Clad.C[0],
                Radi     = FiberA['Clad'][1],
                Index    = FiberA['Clad'][0] )

Clad1 = Circle( Position = Clad.C[0],
                Radi     = FiberA['Clad'][1],
                Index    = FiberA['Clad'][0] )

Clad2 = Circle( Position = Clad.C[0],
                Radi     = FiberB['Clad'][1],
                Index    = FiberB['Clad'][0] )


Core0 = Circle( Position = Clad.C[0],
                Radi     = FiberA['Core'][1],
                Index    = FiberA['Core'][0] )

Core1 = Circle( Position = Clad.C[0],
                Radi     = FiberA['Core'][1],
                Index    = FiberA['Core'][0] )

Core2 = Circle( Position = Clad.C[0],
                Radi     = FiberB['Core'][1],
                Index    = FiberB['Core'][0] )

Geo = Geometry(Objects = [Clad, Clad0, Clad2, Core0, Core1, Core2],
               Xbound  = [-120, 120],
               Ybound  = [-110, 130],
               Nx      = 100,
               Ny      = 100,
               Length  = None)


Sol = SuPySolver(Coupler=Geo)

SuperModes = Sol.GetModes(wavelength = 1.55,
                          Nstep      = 20,
                          Nsol       = 7,
                          debug      = False,
                          ITRi       = 1,
                          ITRf       = 0.05,
                          tolerance  = 1e-20,
                          error      = 3,
                          Xsym       = 0,
                          Ysym       = 0 )

SuperModes[0].PlotPropagation()


SuperModes.SaveFig(Directory  = 'perso/DoubleClad_2',
                   Input      = ['All'],
                   nMax       = 4,
                   iter       = [0,-1])
