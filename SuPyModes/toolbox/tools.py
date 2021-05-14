
#-------------------------Importations------------------------------------------
""" standard imports """
import os
import pandas as pd
import numpy as np

global args
import config.args as args
#-------------------------Importations------------------------------------------


def interactive_print(string, arguments):
        """ This evaluate the 'quiet' status before print are being rendered.

        """
        if args.quiet:
            pass
        else:
            print(string)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def interactive_input(arguments):
    """
    Function that mute inputs to be able to run the unittests
    """

    if args.quiet:
        pass
    else:
        input("Press [enter] to continue.")


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def verify_file(dir):

    if os.path.isfile(dir):
        if args.quiet:
            os.remove(dir)
            return dir

        res = input('File already existing, do you want to rewrite it? (y,n): ')

        if res == 'y':
            os.remove(dir)
        else:
            dir = None

    return dir



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def read_json_file(file_name):

    with open(file_name, 'r') as json_file:

        json_params = json.load(json_file)

        for _, item in enumerate(json_params['geometry']):

            if isinstance(json_params['geometry'][item]['index'], str):
                json_params['geometry'][item]['index'] = eval(json_params['geometry'][item]['index'] )

            json_params['geometry'][item]['function'] = eval(json_params['geometry'][item]['function'] )

        if not json_params['iteration']['variable_list']:
            json_params['iteration']['variable_list'] = np.linspace(json_params['iteration']['start'],
                                                                    json_params['iteration']['stop'],
                                                                    json_params['iteration']['N_step'])

        del json_params['iteration']['start']
        del json_params['iteration']['stop']
        del json_params['iteration']['N_step']

    if json_params['iteration']['symmetric']:
        json_params['iteration']['variable_list'] = np.concatenate((json_params['iteration']['variable_list'][0:-1:1],json_params['iteration']['variable_list'][::-1]))

    return json_params



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def _find_closest_(arr, n, target):

    if (target <= arr[0]):
        return 0
    if (target >= arr[n - 1]):
        return n - 1

    i = 0; j = n; mid = 0
    while (i < j):
        mid = int((i + j) / 2)

        if (arr[mid] == target):
            return mid


        if (target < arr[mid]) :


            if (mid > 0 and target > arr[mid - 1]):
                return getClosest(arr[mid - 1], arr[mid], target)

            j = mid

        else :
            if (mid < n - 1 and target < arr[mid + 1]):
                return getClosest(arr, mid, mid + 1, target)

            i = mid + 1

    return mid


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def _get_closest_(arr, val1, val2, target):

    if (target - arr[val1] >= arr[val2] - target):
        return val2
    else:
        return val1


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


class FiberParameters(object):
    """ This method compute the V value associated to the clad.

    """

    def __init__(self,
                 lamb = None,
                 radius = None,
                 n_clad = None,
                 n_core = None,
                 V_value = None,
                 NA = None):

        self.lamb = lamb
        self.radius = radius
        self.n_clad = n_clad
        self.n_core = n_core
        self.V_value = V_value
        self.fiber = {}
        self.NA = NA


        self.fiber['n_core'] = self.n_core
        self.fiber['NA'] = self.NA
        self.fiber['V_value'] = self.V_value
        self.fiber['n_clad'] = self.n_clad
        self.fiber['radius'] = self.radius
        self.fiber['lambda'] = self.lamb




    def compute_NA(self):
        self.NA = np.sqrt(self.n_core**2 - self.n_clad**2)
        self.fiber['NA'] = self.NA


    def compute_V_value(self):
        self.V_value = 2 * np.pi * self.radius / self.lamb * self.NA
        self.fiber['V_value'] = self.V_value


    def compute_n_clad(self):

        term1 = self.V_value * self.lamb
        term2 = 2 * np.pi * self.radius
        term3 = (term1 / term2)**2
        self.n_clad = ( self.n_core**2 - term3 )**(0.5)
        self.fiber['n_clad'] = self.n_core

    def compute_n_core(self):

        term1 = self.V_value * self.lamb
        term2 = 2*np.pi * self.radius
        term3 = (term1 / term2)**2
        self.n_core = ( self.n_clad**2 + term3 )**(0.5)
        self.fiber['n_core'] = self.n_core







# --
