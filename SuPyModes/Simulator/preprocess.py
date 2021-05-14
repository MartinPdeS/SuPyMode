
#-------------------------Importations------------------------------------------
""" stantdard imports """
import json, os, pickle, sys, itertools
import numpy as np
import pandas as pd

""" package imports """
global args
import SuPyModes.config.args as args
from SuPyModes.toolbox.utils import get_project_root
root = get_project_root()
#-------------------------Importations------------------------------------------

class configuration():

    def __init__(self, recipe_dir):

        with open(args.input_file, 'r') as json_file:

            recipe = json.load(json_file)

        self.recipe = recipe

        self.N_solution =  recipe['iteration']['N_solution']
        self.z_symmetric = recipe['iteration']['symmetric']
        self.N_step =      recipe['iteration']['N_step']

        self.NX_point =    recipe['mesh']['NX_point']
        self.NY_point =    recipe['mesh']['NY_point']
        self.shape =       self.NX_point, self.NY_point
        self.X_bound =     recipe['mesh']["X_boundary"]
        self.Y_bound =     recipe['mesh']["Y_boundary"]
        self.size =        recipe['mesh']['NX_point']*recipe['mesh']['NY_point']

        self.x_symmetry =  recipe['physics']['x_symmetry']
        self.y_symmetry =  recipe['physics']['y_symmetry']
        self.wavelength =  recipe['physics']['lambda']
        self.error_order = recipe['physics']['error order']
        self.tolerance =   recipe['physics']['tolerance']

        if 1 in self.y_symmetry or -1 in self.y_symmetry:
            tempy = self.NY_point * 2
        else:
            tempy = self.NY_point

        if 1 in self.x_symmetry or -1 in self.x_symmetry:
            tempx = self.NX_point * 2
        else:
            tempx = self.NX_point

        self.real_shape = tempx, tempy


        self.degenerate_factor = 1
        self.variable = 'ITR'
        self.recipe_name = recipe_dir

        self.generate_solution_frame()

        self.generate_config_file()


    def generate_solution_frame(self):

        array = [[0,1,-1],[0,1,0],[0,1,1],
                 [0,0,-1],[0,0,0],[0,0,1],
                 [0,-1,-1],[0,-1,0],[0,-1,1]]

        sol_frame = pd.DataFrame( array, columns=['solution','x_symmetry','y_symmetry'])

        sol_frame = sol_frame[sol_frame.x_symmetry.isin( self.recipe['physics']['x_symmetry'] ) ]

        sol_frame = sol_frame[sol_frame.y_symmetry.isin( self.recipe['physics']['y_symmetry'] ) ]


        temp = np.linspace(0,
                           sol_frame.shape[0]-1,
                           self.recipe['iteration']['N_solution']).round().tolist()

        dist = [temp.count(i) for i in set(temp)]

        sol_frame['solution'] = dist

        return sol_frame



    def generate_config_file(self):


        degenerate_factor = 1

        physics = self.recipe['physics']

        mesh = self.recipe['mesh']

        iteration = self.recipe['iteration']

        temp = np.linspace(iteration['start'],iteration['stop'],iteration['N_step'])

        self.variable_list = list(temp)

        frame = self.generate_solution_frame()


        frame.to_pickle(os.path.join(root, 'SuPyModes/temporary_data/solution_frame.pickle'))





        x_vector = np.linspace(self.X_bound, self.NX_point)
        y_vector = np.linspace(self.Y_bound, self.NY_point)
        X_mesh, Y_mesh = np.meshgrid(x_vector, y_vector)
