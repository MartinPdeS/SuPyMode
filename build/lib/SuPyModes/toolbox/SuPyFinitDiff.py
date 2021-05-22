#-------------------------Importations------------------------------------------
""" standard imports """
import numpy as np
from scipy.sparse import spdiags

""" package imports """
global config
import SuPyModes.config.config as config
#-------------------------Importations------------------------------------------

class SuPyFinitdifference(object):


    def __init__(self, Axes):

        self.Axes = Axes


    def second_order_y_symmetry(self, e0, e1, e2, e3, e4):
        """
        Function that return 5 tables, the last of which has been modified
        (part of the coefficients has been multiplied by 2)

        arguments:
            : param e0/e1/e2/e3/e4 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4 : array.

        returns:
            : param e0/e1/e2/e3/e4 : e0/e1/e2/e3 are unchanged and e4 have been modified
            : type e0/e1/e2/e3/e4 : array.

        calls:
            :call1: .second_order_laplacian_sparse
        """

        e4[1:2*self.Axes.Nx] *= 2

        return e0, e1, e2, e3, e4


    def second_order_y_Anti_symmetry(self, e0, e1, e2, e3, e4):
        """
        Function that return 5 tables, the last of which has been modified
        (part of the coefficients has been replaced by 0).

        arguments:
            : param e0/e1/e2/e3/e4 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4 : array.

        returns:
            : param e0/e1/e2/e3/e4 : e0/e1/e2/e3 are unchanged and e4 have been modified
            : type e0/e1/e2/e3/e4 : array.

        calls:
            :call1: XXX
        """
        e4[1:2*self.Axes.Nx] = 0

        return e0, e1, e2, e3, e4


    def third_order_y_symmetry(self, e0, e1, e2, e3, e4, e5, e6, e7, e8):
        """
        Functions that return 9 tables, the last two of which has been modified
        (part of the coefficients has been multiplied by 2)

        arguments:
            : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

        returns:
            : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : e0/e1/e2/e3/e4/e5/e6 are
            unchanged and e7/e8 have been modified
            : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

        calls:
            :call1: .third_order_laplacian_sparse()
        """

        e7[1:2*self.Axes.Nx] *= 2
        e8[1:4*self.Axes.Nx] *= 2

        return e0, e1, e2, e3, e4, e5, e6, e7, e8


    def third_order_y_Anti_symmetry(self, e0, e1, e2, e3, e4, e5, e6, e7, e8):
        """
        Functions that return 9 tables, the last two of which has been modified
        (part of the coefficients has been replaced by 0)

        arguments:
            : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

        returns:
            : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : e0/e1/e2/e3/e4/e5/e6 are
            unchanged and e7/e8 have been modified
            : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

        calls:
            :call1: XXX
        """
        e7[1:2*self.Axes.Nx] = 0
        e8[1:4*self.Axes.Nx] = 0

        return e0, e1, e2, e3, e4, e5, e6, e7, e8



    def second_order_x_symmetry(self, e0, e1, e2, e3, e4):
        """
        Function that return 5 tables, the last of which has been modified
        (part of the coefficients has been multiplied by 2)

        arguments:
            : param e0/e1/e2/e3/e4 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4 : array.

        returns:
            : param e0/e1/e2/e3/e4 : e0/e1/e3/e4 are unchanged and e2 have been modified
            : type e0/e1/e2/e3/e4 : array.

        calls:
            :call1: .second_order_laplacian_sparse
        """
        for j in range(self.Axes.Ny):
            e2[j*self.Axes.Nx+1] *= 2

        return e0, e1, e2, e3, e4


    def second_order_x_Anti_symmetry(self, e0, e1, e2, e3, e4):
        """
        Function that return 5 tables, the last of which has been modified
        (part of the coefficients has been replaced by 0).

        arguments:
            : param e0/e1/e2/e3/e4 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4 : array.

        returns:
            : param e0/e1/e2/e3/e4 : e0/e1/e3/e4 are unchanged and e2 have been modified
            : type e0/e1/e2/e3/e4 : array.

        calls:
            :call1: XXX
        """
        for j in range(self.Axes.Ny):
            e2[j*self.Axes.Nx+1] = 0

        return e0, e1, e2, e3, e4


    def third_order_x_symmetry(self, e0, e1, e2, e3, e4, e5, e6, e7, e8):
        """
        Functions that return 9 tables, the last two of which has been modified
        (part of the coefficients has been multiplied by 2)

            arguments:
                : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : tables from the finite difference coefficients
                : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

            returns:
                : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : e0/e1/e2/e5/e6/e7/e8 are
                unchanged and e3/e4 have been modified
                : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

            calls:
                :call1: .third_order_laplacian_sparse()
        """
        for j in range(self.Axes.Ny):
            e3[j*self.Axes.Nx+1] *= 2
            e4[j*self.Axes.Nx+2] *= 2
            e4[j*self.Axes.Nx+3] *= 2

        return e0, e1, e2, e3, e4, e5, e6, e7, e8


    def third_order_x_Anti_symmetry(self, e0, e1, e2, e3, e4, e5, e6, e7, e8):
        """
        Functions that return 9 tables, the last two of which has been modified
        (part of the coefficients has been replaced by 0)

        arguments:
            : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

        returns:
            : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : e0/e1/e2/e3/e4/e5/e6 are
            unchanged and e7/e8 have been modified
            : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

        calls:
            :call1: XXX

        """
        for j in range(self.Axes.Ny):
            e3[j*self.Axes.Nx+1] = 0
            e4[j*self.Axes.Nx+2] = 0
            e4[j*self.Axes.Nx+3] = 0

        return e0, e1, e2, e3, e4, e5, e6, e7, e8


    def second_order_Dirichlet(self):
        """
        Function that creates 5 tables from the finite difference coefficients
        (second derivative, second order). This method return diags element with
        Dirichlet boundary conditions.

        returns:
            : param e0/e1/e2/e3/e4 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4 : array.

        calls:
            :call1: .second_order_laplacian_sparse()
        """

        e0 = np.array(
                     [2] * self.Axes.Nx * self.Axes.Ny)


        e1 = np.array(
                     ([-1]*(self.Axes.Nx-1) + [0]) * (self.Axes.Ny))


        e2 = np.array(
                     ([0] + [-1]*(self.Axes.Nx-1)) * (self.Axes.Ny))


        e3 = np.array(
                      [-1]*self.Axes.Nx * self.Axes.Ny)


        e4 = np.array(
                     [-1]*self.Axes.Nx * self.Axes.Ny)


        return e0, e1, e2, e3, e4

    def third_order_Dirichlet(self):
        """
        Function that creates 9 tables from the finite difference coefficients
        (second derivative, fourth order).This method return diags element with
        Dirichlet boundary conditions.

        returns:
            : param e0/e1/e2/e3/e4/e5/e6/e7/e8 : tables from the finite difference coefficients
            : type e0/e1/e2/e3/e4/e5/e6/e7/e8 : array.

        calls:
            :call1: .third_order_laplacian_sparse()
        """

        e0 = np.array(
                     [5/2] * self.Axes.Nx * self.Axes.Ny)


        e1 = np.array(
                     ([1/12]*(self.Axes.Nx-2) + 2*[0]) * (self.Axes.Ny))


        e2 = np.array(
                     ([-4/3]*(self.Axes.Nx-1) + [0]) * (self.Axes.Ny))


        e3 = np.array(
                      ([0] + [-4/3]*(self.Axes.Nx-1)) * (self.Axes.Ny))


        e4 = np.array(
                     (2*[0] + [1/12]*(self.Axes.Nx-2)) * (self.Axes.Ny))



        e5 = np.array(
                      [1/12]*self.Axes.Nx * self.Axes.Ny)


        e6 = np.array(
                     [-4/3]*self.Axes.Nx * self.Axes.Ny)


        e7 = np.array(
                      [-4/3]*self.Axes.Nx * self.Axes.Ny)


        e8 = np.array(
                     [1/12]*self.Axes.Nx * self.Axes.Ny)

        return e0, e1, e2, e3, e4, e5, e6, e7, e8


    def laplacian_sparse(self, nk, x_symmetry=False, y_symmetry=False):
        """
        arguments:
            : param nk: Matrix containing the values n**2*k**2 for solving the
            Helmoltz equation
            : type nk: array.

            :param order: Order for finite difference resolution
            :type order: int.

        calls:
            :call1: solver.laplacian_sparse()
        """

        if config.error_order == 2:
            self.second_order_laplacian_sparse(nk, x_symmetry, y_symmetry)

        if config.error_order == 3:
            self.third_order_laplacian_sparse(nk, x_symmetry, y_symmetry)


    def second_order_laplacian_sparse(self, nk, x_symmetry, y_symmetry):
        """
        Construct a sparse matrix that applies the 5-point laplacian discretization

        arguments:

            : param nk: Matrix containing the values n**2*k**2 for solving the
            Helmoltz equation
            : type nk: array.

        calls:
            :call1: .laplacian_sparse()
        """

        e0, e1 , e2, e3, e4 = self.second_order_Dirichlet()

        if y_symmetry == 1:
            e0, e1, e2, e3, e4 = self.second_order_y_symmetry(e0, e1, e2, e3, e4)

        if y_symmetry == -1:
            e0, e1, e2, e3, e4 = self.second_order_y_Anti_symmetry(e0, e1, e2, e3, e4)

        if x_symmetry == 1:
            e0, e1, e2, e3, e4 = self.second_order_x_symmetry(e0, e1, e2, e3, e4)

        if x_symmetry == -1:
            e0, e1, e2, e3, e4 = self.second_order_x_Anti_symmetry(e0, e1, e2, e3, e4)

        self.Axes.Direct.dy, self.Axes.Direct.dx = self.Axes.Direct.dx, self.Axes.Direct.dy

        E0 = e0*(1./self.Axes.Direct.dy**2 + 1./self.Axes.Direct.dx**2) - nk

        E1 = e1/self.Axes.Direct.dx**2
        E2 = e2/self.Axes.Direct.dx**2

        E3 = e3/self.Axes.Direct.dy**2
        E4 = e4/self.Axes.Direct.dy**2

        self.Matrix = spdiags([E1, E2, E0, E3, E4],
                              [-1, 1, 0, -self.Axes.Nx, self.Axes.Nx],
                              self.Axes.Nx*self.Axes.Ny, self.Axes.Nx*self.Axes.Ny)


    def third_order_laplacian_sparse(self, nk, x_symmetry, y_symmetry):
        """
        Construct a sparse matrix that applies the 5-point laplacian discretization.

        arguments:

            : param nk: Matrix containing the values n**2*k**2 for solving the
            Helmoltz equation
            : type nk: array.

        calls:
            :call1: .laplacian_sparse()
        """


        e0, e1, e2, e3, e4, e5, e6, e7, e8 = self.third_order_Dirichlet()

        self.Axes.update_deltas()

        if y_symmetry == 1:
            e0, e1, e2, e3, e4, e5, e6, e7, e8 = self.third_order_y_symmetry(e0, e1, e2, e3, e4, e5, e6, e7, e8)

        if y_symmetry == -1:
            e0, e1, e2, e3, e4, e5, e6, e7, e8 = self.third_order_y_Anti_symmetry(e0, e1, e2, e3, e4, e5, e6, e7, e8)

        if x_symmetry == 1:
            e0, e1, e2, e3, e4, e5, e6, e7, e8 = self.third_order_x_symmetry(e0, e1, e2, e3, e4, e5, e6, e7, e8)

        if x_symmetry == -1:
            e0, e1, e2, e3, e4, e5, e6, e7, e8 = self.third_order_x_Anti_symmetry(e0, e1, e2, e3, e4, e5, e6, e7, e8)


        E0 = e0*(1./self.Axes.Direct.dy**2 + 1./self.Axes.Direct.dx**2) - nk

        E1 = e1/self.Axes.Direct.dx**2
        E2 = e2/self.Axes.Direct.dx**2
        E3 = e3/self.Axes.Direct.dx**2
        E4 = e4/self.Axes.Direct.dx**2

        E5 = e5/self.Axes.Direct.dy**2
        E6 = e6/self.Axes.Direct.dy**2
        E7 = e7/self.Axes.Direct.dy**2
        E8 = e8/self.Axes.Direct.dy**2

        self.Matrix = spdiags([E1, E2, E3, E4, E0, E5, E6, E7, E8],
                              [-2, -1, 1, 2, 0, -2*self.Axes.Nx, -self.Axes.Nx, self.Axes.Nx, 2*self.Axes.Nx],
                              self.Axes.Nx*self.Axes.Ny, self.Axes.Nx*self.Axes.Ny)
