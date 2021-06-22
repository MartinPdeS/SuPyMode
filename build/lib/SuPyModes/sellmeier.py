import numpy as np



def BK7_glass(wavelength):
    B1 = 1.03961212
    B2 = 0.231792344
    B3 = 1.01046945

    C1 = 6.00069867e-3 #[micro M^2]
    C2 = 2.00179144e-2 #[micro M^2]
    C3 = 1.03560653e2 #[micro M^2]

    term0 = 1
    term1 = B1*wavelength**2/(wavelength**2 - C1)
    term2 = B2*wavelength**2/(wavelength**2 - C2)
    term3 = B3*wavelength**2/(wavelength**2 - C3)

    n = np.sqrt(term0 + term1 + term2 + term3)

    return n


def Fused_silica(wavelength):
    B1 = 0.696166300
    B2 = 0.407942600
    B3 = 0.897479400

    C1 = 4.67914826e-3 #[micro M^2]
    C2 = 1.35120631e-2 #[micro M^2]
    C3 = 97.9340025    #[micro M^2]

    term0 = 1
    term1 = B1*wavelength**2/(wavelength**2 - C1)
    term2 = B2*wavelength**2/(wavelength**2 - C2)
    term3 = B3*wavelength**2/(wavelength**2 - C3)

    n = np.sqrt(term0 + term1 + term2 + term3)

    return n


def Ambiant_air(wavelength):
    B1 = 4.915889e-4
    B2 = 5.368273e-5
    B3 = -1.949567e-4

    C1 = 4.352140e-3 #[micro M^2]
    C2 = 1.747001e-2 #[micro M^2]
    C3 = 4.258444e3  #[micro M^2]

    term0 = 1
    term1 = B1*wavelength**2/(wavelength**2 - C1)
    term2 = B2*wavelength**2/(wavelength**2 - C2)
    term3 = B3*wavelength**2/(wavelength**2 - C3)

    n = np.sqrt(term0 + term1 + term2 + term3)

    return n
