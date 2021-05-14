
import pickle
global config
import config.config as config

def Step_func(**kwargs):
    """
    This function compute a 2D smooth symmetrique step function.

    arguments:
        : param x: x value in the mesh
        : type x: float.

        : param y: y value in the mesh
        : type y: float.

        : param radius: the radius of the step function
        : type radius: float.

        : param x_center: the center in the x axis of the step function
        : type x_center: float.

        : param y_center: the center in the y axis of the step function
        : type y_center: float.

    returns:
        : param val:
        : type val: bool.

    calls:
        :call1:

    """

    with open( "src/python/config/mesh.pickle", "rb" ) as f:
        mesh = pickle.load(f)

    return (((mesh['X'].T-kwargs['x_center'])**2 + (mesh['Y'].T-kwargs['y_center'])**2) <= kwargs['radius']**2)



def Step_modified_func(**kwargs):
    """
    This function compute a 2D smooth symmetrique step function.

    arguments:
        : param x: x value in the mesh
        : type x: float.

        : param y: y value in the mesh
        : type y: float.

        : param radius: the radius of the step function
        : type radius: float.

        : param x_center: the center in the x axis of the step function
        : type x_center: float.

        : param y_center: the center in the y axis of the step function
        : type y_center: float.

    returns:
        : param val:
        : type val: bool.

    calls:
        :call1:

    """

    with open( "src/python/config/mesh.pickle", "rb" ) as f:
        mesh = pickle.load(f)

    return ((((mesh['X']-kwargs['x_center'])/0.95)**2 + (mesh['Y']-kwargs['y_center'])**2) <= kwargs['radius']**2)


def Load_profil(dir):
    """
    This method load a new image as index profil.

    arguments:
    :param dir: Directory of image to load (default=uint8).
    :type dir: str.

    :param Nx: Length of the output matrix.
    :type Nx: int.

    :param Ny: Width of the output matrix.
    :type Ny: int.

    returns:
        : param image:
        : type image: array.

    calls:
        :call1:

    """

    dir = index_dir + dir

    image = cv2.imread(dir, 0)

    image = cv2.resize(image,
                      (config.Ny, config.Nx),
                      interpolation = cv2.INTER_LINEAR)

    return image



# --
