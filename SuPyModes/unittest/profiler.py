
#-------------------------Importations------------------------------------------
""" stantdard imports """
import os, sys


""" package imports """
from SuPyModes.Simulator.solver import SuPySolver
from SuPyModes.Simulator.post_processing import SuPyProcessing
from SuPyModes.Simulator.data_collector import collector
from SuPyModes.toolbox.utils import get_project_root
from SuPyModes.Simulator.geometry import Coupler2
from SuPyModes.Simulator.preprocess import configuration
from SuPyModes.toolbox.utils import get_project_root
root = get_project_root()
from SuPyModes.toolbox.argument_parser import parse
from SuPyModes.toolbox.sellmeier import Fused_silica
#-------------------------Importations------------------------------------------

args = parse()

def main():

    config = configuration(args.input_file)

    Obj = Coupler2()

    Obj.add_clad(init=200,
                 R_clad0=62.5,
                 R_clad1=62.5,
                 fusion=0.9,
                 index=Fused_silica(1.55))

    Obj.add_cores(position='core0',
                  radius=4.1,
                  index=Fused_silica(1.55)+0.005)

    Obj.add_cores(position='core1',
                  radius=4.1,
                  index=Fused_silica(1.55)+0.005)

    Obj.rasterize_mesh([0,100], [0,100], 200,200)

    #Obj.plot_geometry()

    Sol = SuPySolver(coupler=Obj)

    Sol.main_loop(wavelength=1.55,
                  Nstep=20,
                  Nsol=5,
                  Xsym=1,
                  Ysym=1)

    Post = SuPyProcessing(Obj.mesh, Sol.Fields)

    Post.compute_coupling()

    Coll = collector(Post.Fields)

    Coll.compute()


if __name__ == "__main__":
    if args.profiler:

        from pycallgraph import PyCallGraph
        from pycallgraph.output import GraphvizOutput
        with PyCallGraph(output=GraphvizOutput()):
            main()
        sys.stdout.write('Profiler complete\n')



    else:
        main()



# -
