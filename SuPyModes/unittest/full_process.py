
#-------------------------Importations------------------------------------------
""" stantdard imports """
import os, sys


""" package imports """
from SuPyModes.Simulator.solver import SuPySolver
from SuPyModes.Simulator.post_processing import SuPyProcessing
from SuPyModes.Simulator.data_collector import collector
from SuPyModes.Simulator.geometry import Coupler2

from SuPyModes.toolbox.sellmeier import Fused_silica
from SuPyModes.toolbox.utils import get_project_root
root = get_project_root()

#-------------------------Importations------------------------------------------

def main():

    Coupler = Coupler2()

    Coupler.add_clad(init=200,
                 R_clad0=62.5,
                 R_clad1=62.5,
                 fusion=0.99,
                 index=Fused_silica(1.55))

    Coupler.add_cores(position='core0',
                  radius=4.1,
                  index=Fused_silica(1.55)+0.005
                  )

    Coupler.add_cores(position='core1',
                  radius=4.1,
                  index=Fused_silica(1.55)+0.005
                  )

    Coupler.rasterize_mesh(Xbound=[0,100],
                           Ybound=[0,100],
                           Nx=60,
                           Ny=60)

    Data = SuPySolver(coupler=Coupler)

    Data.main_loop(wavelength=1.55,
                  Nstep=20,
                  ITRf=0.05,
                  debug=False,
                  naming=True,
                  Nsol=7,
                  Xsym=1,
                  Ysym=1)

    Post = SuPyProcessing(Data)

    Post.compute_coupling()

    Coll = collector(Post.Fields)

    Coll.compute()


if __name__ == "__main__":

        main()




# -
