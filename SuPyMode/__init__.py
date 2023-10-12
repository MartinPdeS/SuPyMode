import pickle
from pathlib import Path
from SuPyMode.tools.directories import instance_directory


def load_superset(filename: str, directory: str = '.'):
    """
    Saves the superset instance as a serialized pickle file.

    :param      filename:  The filename
    :type       filename:  str
    """
    if directory == 'auto':
        directory = instance_directory

    filename = Path(directory).joinpath(filename).with_suffix('.pickle')

    with open(filename, 'rb') as input_file:
        superset = pickle.load(input_file)

    return superset

# -
