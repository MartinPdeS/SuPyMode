from SuPyModes.Geometry          import Coupler2, Coupler1
from SuPyModes.Solver            import SuPySolver
from SuPyModes.toolbox.sellmeier import Fused_silica

Coupler = Coupler2()


Coupler.add_clad(init=200, R_clad0=62.5, R_clad1=62.5, fusion=0.99, index=Fused_silica(1.55) )

Coupler.add_cores(position='core0', radius=4.1, index=Fused_silica(1.55)+0.005 )

Coupler.add_cores(position='core1', radius=4.1, index=Fused_silica(1.55) + 0.005 )

Coupler.CreateMesh(Xbound=[0,100], Ybound=[0,100], Nx=60, Ny=60)


#Coupler.Plot()

Sol = SuPySolver(coupler=Coupler)

SuperModes = Sol.GetModes(wavelength=1.55,
                          Nstep=100,
                          Nsol=3,
                          Xsym=1,
                          Ysym=1)

SuperModes.Plot('Coupling',iter=0)
