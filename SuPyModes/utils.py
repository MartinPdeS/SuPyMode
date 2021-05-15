import numpy as np


def RecomposeSymmetries(Input, Symmetries, Xaxis=None, Yaxis=None):

    if Symmetries[1] == 1:
        Input = np.concatenate((Input[::-1,:],Input),axis=0)

    if Symmetries[1] == -1:
        Input = np.concatenate((-Input[::-1,:],Input),axis=0)

    if Symmetries[0] == 1:
        Input = np.concatenate((Input[:,::-1],Input),axis=1)

    if Symmetries[0] == -1:
        Input = np.concatenate((-Input[:,::-1],Input),axis=1)

    if Xaxis is not None and Symmetries[0] != 0:
        Xaxis = np.concatenate( (-Xaxis[::-1],Xaxis) )

    if Yaxis is not None and Symmetries[1] != 0:
        Yaxis = np.concatenate( (-Yaxis[::-1],Yaxis) )

    return Input, Xaxis, Yaxis


def CheckSymmetries(SuperMode0, SuperMode1):
    if SuperMode0.xSym[0] == 0 or SuperMode1.xSym[0] == 0: return True

    if SuperMode0.ySym[0] == 0 or SuperMode1.ySym[0] == 0: return True

    if SuperMode0.xSym[0] == - SuperMode1.xSym[0]: return False

    if SuperMode0.ySym[0] == - SuperMode1.ySym[0]: return False

    return True
