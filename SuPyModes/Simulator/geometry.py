#-------------------------Importations------------------------------------------
""" standard imports """
import pickle, json, sys
import numpy as np
import os
from matplotlib.path import Path
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.optimize import least_squares
pi, cos, sin, sqrt = np.pi, np.cos, np.sin, np.sqrt

import shapely.geometry as geometry
from shapely.geometry import Polygon
from shapely.geometry.point import Point
import shapely.ops as so
from shapely.geometry import LineString
from shapely.ops import cascaded_union
from descartes import PolygonPatch
from shapely.ops import nearest_points

""" package imports"""
from SuPyModes.toolbox.utils import get_project_root
root = get_project_root()

#-------------------------Importations------------------------------------------


class Coupler2(object):

    def __init__(self, debug: bool = False):

        self.debug = debug
        self.GEO = None
        self.OBJECT = {}
        self.N = 0
        self.cst = 2*(2-sqrt(2))



    def plot_geometry(self):

        fig = plt.figure(figsize=(15,5))

        ax = fig.add_subplot(131)
        ax1 = fig.add_subplot(132)
        ax2 = fig.add_subplot(133)

        ax.grid()
        ax1.grid()

        ax.set_ylabel(r'Distance Y-axis [$\mu m$]')
        ax.set_xlabel(r'Distance X-axis [$\mu m$]')
        ax.set_title(r'Geometrical configuration', fontsize=10)


        [ ax.plot(*object['shape'].exterior.xy) for object in self.OBJECT.values()]

        [ ax.plot(*object.exterior.xy) for object in self.GEO]

        [ ax.scatter(C.x, C.y) for C in self.C]
        [ ax.text(C.x,C.y,'C{0}'.format(n+1)) for n, C in enumerate(self.C) ]


        [ ax.scatter(P.x, P.y) for P in self.P]
        [ ax.text(P.x, P.y,'P{0}'.format(n+1)) for n, P in enumerate(self.P) ]

        [ ax.scatter(F.x, F.y) for F in self.F]
        [ ax.text(F.x, F.y,'F{0}'.format(n+1)) for n, F in enumerate(self.F) ]




        if not self.coupler.is_empty:

            ax1.add_patch(PolygonPatch(self.coupler, alpha=0.5,))
            ax1.plot()
            ax1.grid()
            ax1.set_xlim([-2*self.R_clad[0],2*self.R_clad[0]])
            ax1.set_ylim([-2*self.R_clad[0],2*self.R_clad[0]])

            ax1.set_title('Fusion degree:{0:.3f}'.format(self.f), fontsize=10)





        [ax1.add_patch(PolygonPatch(object['shape'], alpha=0.5)) for object in self.OBJECT.values() ]



        ax2.set_title('Rasterized RI profil', fontsize=10)

        pcm = ax2.pcolormesh(np.linspace(*self.Ybound, np.shape(self.mesh)[1]),
                             np.linspace(*self.Xbound, np.shape(self.mesh)[0]),
                             self.mesh,
                             cmap='PuBu_r',
                             norm=colors.LogNorm(vmin = 1.42, vmax=self.mesh.max()))

        plt.show()


    def compute_limit(self):

        temp =  so.cascaded_union([self.Clad[0], self.Clad[1]])

        return temp.convex_hull.difference(temp).area


    def construct_fibers(self):

        self.F, self.Clad = [], []

        self.F = [Point(-self.delta/2, 0), Point(self.delta/2, 0)]

        self.Clad = ( self.F[0].buffer(self.R_clad[0]), \
                      self.F[1].buffer(self.R_clad[1]) )

        self.F += ( self.Clad[0].exterior.intersection(self.Clad[1].exterior)[0], \
                    self.Clad[0].exterior.intersection(self.Clad[1].exterior)[1] )


        self.overlap = self.Clad[0].intersection(self.Clad[1])

        self.topo = 'concave' if self.overlap.area > self.compute_limit() else 'convex'



    def find_virtual_coord_from_Rv(self, Rv):


        if self.topo == 'convex':
            temp0 = Point(self.F[0]).buffer(self.R_clad[0] + Rv)

            temp1 = Point(self.F[1]).buffer(self.R_clad[1] + Rv)

        else:

            temp0 = Point(self.F[0]).buffer(Rv - self.R_clad[0])

            temp1 = Point(self.F[1]).buffer(Rv - self.R_clad[1])


        self.C += ( temp0.exterior.intersection(temp1.exterior)[0], \
                    temp0.exterior.intersection(temp1.exterior)[1] )



    def construct_triangle_from_Rv(self, Rv):

        if self.topo == 'concave' and Rv < self.R_clad[0] :

            Rv = self.R_clad[0]*1.1

        self.find_virtual_coord_from_Rv(Rv)

        self.triangle_top, self.triangle_bottom = Polygon([self.F[0], self.F[1], self.C[0]]), Polygon([self.F[0], self.F[1], self.C[1]])

        self.GEO += [self.triangle_top, self.triangle_bottom]

        return Rv




    def construct_virtual(self, Rv):

        self.circ2, self.circ3  = self.C[0].buffer(Rv), self.C[1].buffer(Rv)

        self.GEO += [self.circ2, self.circ3]

        self.P = self.P + [nearest_points(self.Clad[1].exterior, self.circ3.exterior)[1],
                           nearest_points(self.Clad[0].exterior, self.circ3.exterior)[1],
                           nearest_points(self.Clad[0].exterior, self.circ2.exterior)[1],
                           nearest_points(self.Clad[1].exterior, self.circ2.exterior)[1],
                           ]



    def compute_added(self):

        if self.topo == 'convex':

            union_top = so.cascaded_union([self.Clad[0],self.Clad[1], self.circ2])

            union_bottom = so.cascaded_union([self.Clad[0],self.Clad[1], self.circ3])

            added_top = self.triangle_top.difference(union_top)

            added_bottom = self.triangle_bottom.difference(union_bottom)

            self.added = so.cascaded_union([added_top, added_bottom])

            self.coupler = so.cascaded_union([self.added, self.Clad[0], self.Clad[1]])


        else:

            Circles = so.cascaded_union([ self.Clad[0], self.Clad[1]])

            envelop = self.circ2.intersection(self.circ3).intersection(self.mask)

            self.coupler = so.cascaded_union([ Circles, envelop])

            self.added = self.coupler.difference(Circles)



    def create_mask(self):

        P1 = nearest_points(self.Clad[1].exterior,self.circ2.exterior)

        P2 = nearest_points(self.Clad[0],self.circ2.exterior)

        P3 = nearest_points(self.Clad[1].exterior,self.circ3.exterior)

        P4 = nearest_points(self.Clad[0].exterior,self.circ3.exterior)

        P1 = Point(P1[1].x, 1000)
        P2 = Point(P2[1].x, 1000)
        P3 = Point(P3[1].x, -1000)
        P4 = Point(P4[1].x, -1000)

        pointList = [P1, P2, P3, P4]

        self.mask = Polygon(pointList).convex_hull



    def fusion_cost(self, Rv):

        self.C, self.P = [], []

        self.GEO = [self.Clad[0], self.Clad[1]]

        Rv = self.construct_triangle_from_Rv(Rv)

        self.construct_virtual(Rv)

        self.create_mask()

        self.compute_added()

        cost = np.abs( self.added.area - self.overlap.area )

        if self.debug:
            sys.stdout.write('Topology:{0} \t Cost:{1:.5f} \t Rv:{2:.0f}\n'.format(self.topo, cost, float(Rv)))

        return cost




    def add_cores(self, position, radius, index=1):

        if position == 'center':
            position = (0,0)

        elif position == 'core0':
            position = self.core_position[0]

        elif position == 'core1':
            position = self.core_position[1]


        shape = Point(position).buffer(radius)

        self.OBJECT[self.N] = {'shape': shape,
                                'index': index}




    def measure_side_area(self, centers):

        p1_right = Point(centers[1], self.Ybound[0])

        p2_right = Point(centers[1], self.Ybound[1])

        p3_right = Point(self.Xbound[1], self.Ybound[1])

        p4_right = Point(self.Xbound[1], self.Ybound[0])


        p1_left = Point(centers[0], self.Ybound[0])

        p2_left = Point(centers[0], self.Ybound[1])

        p3_left = Point(self.Xbound[0], self.Ybound[1])

        p4_left = Point(self.Xbound[0], self.Ybound[0])


        left_point_list = [p1_left, p2_left, p3_left, p4_left]

        right_point_list = [p1_right, p2_right, p3_right, p4_right]

        left_mask = Polygon([[p.x, p.y] for p in left_point_list])

        right_mask = Polygon([[p.x, p.y] for p in right_point_list])

        cost0 = np.abs(self.coupler.intersection(left_mask).area - pi *self.R_clad[0]**2/2)

        cost1 = np.abs(self.coupler.intersection(right_mask).area - pi *self.R_clad[1]**2/2)

        if self.debug:
            sys.stdout.write('Cost:{0:.2f} \t position left:{1:.2f} \t position right:{2:.2f}\n'.format(cost0 + cost1, centers[0], centers[1]))

        return cost1, cost0


    def add_clad(self, init, R_clad0, R_clad1, fusion=0.9, index=1):

        self.f = fusion

        self.R_clad = R_clad0, R_clad1

        self.Xbound = [-2*max(self.R_clad[0],self.R_clad[1]),2*max(self.R_clad[0],self.R_clad[1])]

        self.Ybound = [-max(self.R_clad[0],self.R_clad[1]),max(self.R_clad[0],self.R_clad[1])]

        self.delta = self.R_clad[0] + self.R_clad[1] - 2 * self.f * (self.R_clad[0] + self.R_clad[1] - np.sqrt(self.R_clad[0]**2 + self.R_clad[1]**2) )

        self.construct_fibers()

        if self.topo == 'concave':
            bounds = ((self.R_clad[0]+self.delta+self.R_clad[1])/2,10000)
        else:
            bounds = (0,10000)

        z = least_squares(self.fusion_cost,init,bounds=bounds)

        self.OBJECT['clad'] = {'shape': self.coupler,
                                "index": index}

        self.optimize_core_position()


    def optimize_core_position(self, init=0):


        z = fsolve(self.measure_side_area, [self.F[0].x, self.F[1].x])

        self.core_position = list( zip( z,[0,0] ) )



    def add_object_to_mesh(self, mesh, obj, index):

        mesh *= (-1 * obj + 1)

        mesh += obj * index

        return mesh



    def rasterize_polygone(self, polygone, points, Nx, Ny):

        obj_ext = Path(list( polygone.exterior.coords))

        obj_ext = obj_ext.contains_points(points).reshape([Nx, Ny])

        return obj_ext



    def rasterize_mesh(self, Xbound, Ybound, Nx, Ny):

        self.mesh = np.ones((Nx, Ny))

        self.Xbound, self.Ybound = Xbound, Ybound

        xv = np.linspace(*Xbound,Nx)
        yv = np.linspace(*Ybound,Ny)

        x, y = np.meshgrid(xv,yv)

        x, y = x.flatten(), y.flatten()

        points = np.vstack((x,y)).T

        for k, v in self.OBJECT.items():

            if self.debug:
                sys.stdout.write('Processing:{0} \t'.format(k))

            obj = self.rasterize_polygone(v['shape'], points, Nx, Ny)

            self.mesh = self.add_object_to_mesh(self.mesh, obj, v['index'])

        self.info = {
                     'Xbound': Xbound,
                     'Ybound': Ybound,
                     'Nx': Nx,
                     'Ny': Ny,
                     'shape': [Nx, Ny],
                     'size': np.size(self.mesh),
                     'x_symmetry': None,
                     'y_symmetry': None
                     }



    def save_data(self, dir: str):
        """
        This method update the class Data with new data_set.

        arguments:
            :param dir: Directory to save geometry file.
            :type dir: str.

        calls:
            :call1: processing.process_geometry()

        """

        with open(dir, 'wb+') as f:

            pickle.dump(self.mesh.T, f)






class Coupler1(object):

    def __init__(self, debug: bool = False):

        self.OBJECT = {}

        self.N = 0

        self.debug = debug



    def plot_geometry(self):

        fig = plt.figure(figsize=(15,5))

        ax = fig.add_subplot(131)
        ax1 = fig.add_subplot(132)
        ax2 = fig.add_subplot(133)

        ax.grid()
        ax1.grid()

        ax.set_ylabel(r'Distance Y-axis [$\mu m$]')
        ax.set_xlabel(r'Distance X-axis [$\mu m$]')
        ax.set_title(r'Geometrical configuration', fontsize=10)


        [ ax.plot(*object['shape'].exterior.xy) for object in self.OBJECT.values()]


        [ax1.add_patch(PolygonPatch(object['shape'], alpha=0.5)) for object in self.OBJECT.values() ]




        ax1.plot()
        ax1.grid()
        ax1.set_xlim([-2*self.R_clad,2*self.R_clad])
        ax1.set_ylim([-2*self.R_clad,2*self.R_clad])



        ax2.set_title('Rasterized RI profil', fontsize=10)

        pcm = ax2.pcolormesh(np.linspace(*self.Ybound, np.shape(self.mesh)[1]),
                             np.linspace(*self.Xbound, np.shape(self.mesh)[0]),
                             self.mesh,
                             cmap='PuBu_r',
                             norm=colors.LogNorm(vmin = 1.42, vmax=self.mesh.max()))

        plt.show()





    def add_cores(self, position, radius, index=1):

        if position == 'center':
            position = (0,0)


        shape = Point(position).buffer(radius)

        self.OBJECT[self.N] = {'shape': shape,
                                'index': index}



    def add_clad(self, R_clad, index=1):

        shape = Point([0,0]).buffer(R_clad)
        self.OBJECT['clad'] = {'shape': shape,
                                'index': index}
        self.R_clad = R_clad


    def optimize_core_position(self, init=0):


        z = fsolve(self.measure_side_area, [self.F[0].x, self.F[1].x])

        self.core_position = list( zip( z,[0,0] ) )



    def add_object_to_mesh(self, mesh, obj, index):

        mesh *= (-1 * obj + 1)

        mesh += obj * index

        return mesh



    def rasterize_polygone(self, polygone, points, Nx, Ny):

        obj_ext = Path(list( polygone.exterior.coords))

        obj_ext = obj_ext.contains_points(points).reshape([Nx, Ny])

        return obj_ext



    def rasterize_mesh(self, Xbound, Ybound, Nx, Ny):

        self.Xbound, self.Ybound = Xbound, Ybound

        self.mesh = np.ones((Nx, Ny))

        xv = np.linspace(*Xbound,Nx)
        yv = np.linspace(*Ybound,Ny)

        x, y = np.meshgrid(xv,yv)

        x, y = x.flatten(), y.flatten()

        points = np.vstack((x,y)).T

        for k, v in self.OBJECT.items():

            if self.debug:
                sys.stdout.write('Processing:{0} \t'.format(k))

            obj = self.rasterize_polygone(v['shape'], points, Nx, Ny)

            self.mesh = self.add_object_to_mesh(self.mesh, obj, v['index'])

        self.info = {
                     'Xbound': Xbound,
                     'Ybound': Ybound,
                     'Nx': Nx,
                     'Ny': Ny,
                     'shape': [Nx, Ny],
                     'size': np.size(self.mesh),
                     'x_symmetry': None,
                     'y_symmetry': None
                     }



    def save_data(self, dir: str):
        """
        This method update the class Data with new data_set.

        arguments:
            :param dir: Directory to save geometry file.
            :type dir: str.

        calls:
            :call1: processing.process_geometry()

        """

        with open(dir, 'wb+') as f:

            pickle.dump(self.mesh.T, f)
