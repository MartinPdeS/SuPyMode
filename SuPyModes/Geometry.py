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
from SuPyModes.Directories import RootPath



#-------------------------Importations------------------------------------------

def Overlap(object0, object1):
    return object0.intersection(object1)

def Intersection(object0, object1):
    return object0.exterior.intersection(object1.exterior)

class Coupler2(object):

    def __init__(self, kwargs = {}):

        self.debug = False
        self.N = 0
        self.cst = 2*(2-sqrt(2))
        self.Points    = {'P': [], 'F': [], 'C': []}
        self.Triangles = {}
        self.Circles   = {}
        self.Cladding  = {}
        self.Cores     = {}
        self.Geometry  = {}
        self.InitClad  = {}
        self.BuildGeo  = {}



    def Plot(self):

        fig = plt.figure(figsize=(15,5))

        ax = fig.add_subplot(131)
        ax1 = fig.add_subplot(132)
        ax2 = fig.add_subplot(133)

        ax.grid()
        ax1.grid()

        ax.set_ylabel(r'Distance Y-axis [$\mu m$]')
        ax.set_xlabel(r'Distance X-axis [$\mu m$]')
        ax.set_title(r'Geometrical configuration', fontsize=10)

        [ax.plot(*object.exterior.xy) for object in self.Triangles.values()]

        [ax.plot(*object.exterior.xy) for object in self.Circles.values()]

        for key, object in self.Geometry.items():
            ax1.add_patch(PolygonPatch(object['shape'], alpha=0.5))


        for key, points in self.Points.items():
            [ ax.scatter(p.x, p.y) for p in points]
            [ ax.text(p.x,p.y,f'{key}{n+1}') for n, p in enumerate(points) ]


        if not self.coupler.is_empty:

            ax1.add_patch(PolygonPatch(self.coupler, alpha=0.5,))
            ax1.plot()
            ax1.grid()
            ax1.set_xlim([-2*self.R_clad[0],2*self.R_clad[0]])
            ax1.set_ylim([-2*self.R_clad[0],2*self.R_clad[0]])
            ax1.set_title('Fusion degree:{0:.3f}'.format(self.f), fontsize=10)

        ax2.set_title('Rasterized RI profil', fontsize=10)

        pcm = ax2.pcolormesh(np.linspace(*self.Ybound, np.shape(self.mesh)[1]),
                             np.linspace(*self.Xbound, np.shape(self.mesh)[0]),
                             self.mesh,
                             cmap='PuBu_r',
                             shading='auto',
                             norm=colors.LogNorm(vmin = 1.42, vmax=self.mesh.max()))

        plt.show()


    def compute_limit(self):

        temp =  so.cascaded_union([self.Clad[0], self.Clad[1]])

        return temp.convex_hull.difference(temp).area


    def construct_fibers(self):

        self.Clad = []

        F = [Point(-self.delta/2, 0), Point(self.delta/2, 0)]

        self.Clad = ( F[0].buffer(self.R_clad[0]),
                      F[1].buffer(self.R_clad[1]) )

        self.BuildGeo['clad'] = ( F[0].buffer(self.R_clad[0]),
                                  F[1].buffer(self.R_clad[1]) )

        self.BuildGeo['overlap'] = Overlap(self.BuildGeo['clad'][0],
                                           self.BuildGeo['clad'][1])

        self.Points['F'] = ( *F,
                             *Intersection( self.BuildGeo['clad'][0],
                                            self.BuildGeo['clad'][1] ) )

        self.topo = 'concave' if self.BuildGeo['overlap'].area > self.compute_limit() else 'convex'



    def find_virtual_coord_from_Rv(self, Rv):

        if self.topo == 'convex':
            p0 = Point(self.Points['F'][0]).buffer(self.R_clad[0] + Rv)

            p1 = Point(self.Points['F'][1]).buffer(self.R_clad[1] + Rv)

        else:

            p0 = Point(self.Points['F'][0]).buffer(Rv - self.R_clad[0])

            p1 = Point(self.Points['F'][1]).buffer(Rv - self.R_clad[1])

        self.Points['C'] = Intersection(p0, p1)


    def construct_triangle_from_Rv(self, Rv):

        if self.topo == 'concave' and Rv < self.R_clad[0] :
            Rv = self.R_clad[0]*1.1

        self.find_virtual_coord_from_Rv(Rv)

        self.Triangles['top']    = Polygon([self.Points['F'][0],
                                            self.Points['F'][1],
                                            self.Points['C'][0]])

        self.Triangles['bottom'] = Polygon([self.Points['F'][0],
                                            self.Points['F'][1],
                                            self.Points['C'][1]])

        return Rv




    def construct_virtual(self, Rv):

        self.Circles['2'] = self.Points['C'][0].buffer(Rv)

        self.Circles['3'] = self.Points['C'][1].buffer(Rv)

        self.Points['P']  = [nearest_points(self.Clad[1].exterior, self.Circles['3'].exterior)[1],
                             nearest_points(self.Clad[0].exterior, self.Circles['3'].exterior)[1],
                             nearest_points(self.Clad[0].exterior, self.Circles['2'].exterior)[1],
                             nearest_points(self.Clad[1].exterior, self.Circles['2'].exterior)[1] ]


    def compute_added(self):

        if self.topo == 'convex':
            union_top = so.cascaded_union([self.Clad[0],self.Clad[1], self.Circles['2']])

            union_bottom = so.cascaded_union([self.Clad[0],self.Clad[1], self.Circles['3']])

            added_top = self.self.Triangles['top'].difference(union_top)

            added_bottom = self.self.Triangles['bottom'].difference(union_bottom)

            self.added = so.cascaded_union([added_top, added_bottom])

            self.coupler = so.cascaded_union([self.added, self.Clad[0], self.Clad[1]])


        else:
            Circles = so.cascaded_union([ self.Clad[0], self.Clad[1]])

            envelop = self.Circles['2'].intersection(self.Circles['3']).intersection(self.mask)

            self.coupler = so.cascaded_union([ Circles, envelop])

            self.added = self.coupler.difference(Circles)



    def create_mask(self):

        P1 = nearest_points(self.Clad[1].exterior,self.Circles['2'].exterior)

        P2 = nearest_points(self.Clad[0],self.Circles['2'].exterior)

        P3 = nearest_points(self.Clad[1].exterior,self.Circles['3'].exterior)

        P4 = nearest_points(self.Clad[0].exterior,self.Circles['3'].exterior)

        P1 = Point(P1[1].x, 1000)
        P2 = Point(P2[1].x, 1000)
        P3 = Point(P3[1].x, -1000)
        P4 = Point(P4[1].x, -1000)

        pointList = [P1, P2, P3, P4]

        self.mask = Polygon(pointList).convex_hull



    def fusion_cost(self, Rv):

        self.InitClad['0'] = {'shape': self.Clad[0], 'index': None}

        self.InitClad['1'] = {'shape': self.Clad[1], 'index': None}

        Rv = self.construct_triangle_from_Rv(Rv)

        self.construct_virtual(Rv)

        self.create_mask()

        self.compute_added()

        cost = np.abs( self.added.area - self.BuildGeo['overlap'].area )

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

        self.Cores[self.N] = {'shape': shape, 'index': index}

        self.Geometry[f'Core {self.N}'] = {'shape': shape, 'index': index}


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

        self.Cladding = {'shape': self.coupler, "index": index}

        self.Geometry['clad'] = {'shape': self.coupler, "index": index}

        self.optimize_core_position()


    def optimize_core_position(self, init=0):

        z = fsolve(self.measure_side_area, [self.Points['F'][0].x, self.Points['F'][1].x])

        self.core_position = list( zip( z,[0,0] ) )


    def add_object_to_mesh(self, mesh, obj, index):

        mesh *= (-1 * obj + 1)

        mesh += obj * index

        return mesh


    def rasterize_polygone(self, polygone, points, Nx, Ny):

        obj_ext = Path(list( polygone.exterior.coords))

        obj_ext = obj_ext.contains_points(points).reshape([Nx, Ny])

        return obj_ext


    def CreateMesh(self, Xbound, Ybound, Nx, Ny):

        self.mesh = np.ones((Nx, Ny))

        self.Xbound, self.Ybound = Xbound, Ybound

        xv = np.linspace(*Xbound,Nx)
        yv = np.linspace(*Ybound,Ny)

        x, y = np.meshgrid(xv,yv)

        x, y = x.flatten(), y.flatten()

        points = np.vstack((x,y)).T

        for k, v in self.Geometry.items():
            obj = self.rasterize_polygone(v['shape'], points, Nx, Ny)
            self.add_object_to_mesh(self.mesh, obj, v['index'])


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



    def Plot(self):

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
        self.Cladding = {'shape': shape, 'index': index}
        self.R_clad = R_clad


    def optimize_core_position(self, init=0):


        z = fsolve(self.measure_side_area, [self.Points['F'][0].x, self.Points['F'][1].x])

        self.core_position = list( zip( z,[0,0] ) )



    def add_object_to_mesh(self, mesh, obj, index):

        self.mesh *= (-1 * obj + 1)

        self.mesh += obj * index



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

        for k, v in self.Geometry.items():

            if self.debug:
                sys.stdout.write('Processing:{0} \t'.format(k))

            obj = self.rasterize_polygone(v['shape'], points, Nx, Ny)

            self.add_object_to_mesh(self.mesh, obj, v['index'])

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
