
import sys
import argparse

def parse():
    
    parser = argparse.ArgumentParser(description='Fiber Mode Simulator')

    parser.add_argument('--input-json-file',
                        help='Input json file',
                        dest='input_file',
                        type=str,
                        required=False)

    parser.add_argument('-q', '--quiet', help='cut off verbose',
                        dest='quiet',
                        action='store_true',
                        required=False)

    parser.add_argument('-n', '--name', help='Naming of the modes',
                        dest='name',
                        action='store_true',
                        required=False)

    parser.add_argument('-p', '--plot', help='plot mode at each iteration of main',
                        dest='plot',
                        action='store_true',
                        required=False)

    parser.add_argument('-s', '--save-plots', help='save plots',
                        dest='savefig',
                        action='store_true',
                        required=False)

    parser.add_argument('-d', '--debug', help='debug mode',
                        dest='debug',
                        action='store_true',
                        required=False)

    parser.add_argument('--profile', help='profiler mode',
                        dest='profiler',
                        action='store_true',
                        required=False)


    parser.add_argument('-l', '--coupler-length', help='coupler length in micro-meter',
                        dest='length',
                        type=float,
                        required=False)

    return parser.parse_args()

"""
with open('SuPyModes/config/args.py', 'w+') as f:

    f.write("input_file='{0}'\nquiet={1}\nplot={2}\nsavefig={3}\ndebug={4}\nprofiler={5}\nlength={6}\nname={7}\n"\
    .format(arg.input_file,\
            arg.quiet,\
            arg.plot,\
            arg.savefig,\
            arg.debug,\
            arg.profiler,\
            arg.length,
            arg.name))"""
