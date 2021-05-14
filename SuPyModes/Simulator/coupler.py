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
from SuPyModes.toolbox.sellmeier import *
from SuPyModes.toolbox.utils import get_project_root
root = get_project_root()

global config, args
#import SuPyModes.config.config as config
import SuPyModes.config.args as args
#-------------------------Importations------------------------------------------


class OptGeo(object):

    def __init__(self, config):

        self.OBJECT = None

        self.dict = config.recipe

        self.f = self.dict['geometry']['Clad']["degree"]

        self.N_fiber = self.dict['geometry']['Clad']["fusion"]


        if self.N_fiber == 1:
            self.R_clad = ( self.dict['geometry']['Clad']["radius0"], )

            self.Xbound = [-2*self.R_clad[0],2*self.R_clad[0]]

            self.Ybound = [-self.R_clad[0],self.R_clad[0]]


        if self.N_fiber == 2:
            self.R_clad = (self.dict['geometry']['Clad']["radius0"], self.dict['geometry']['Clad']["radius1"])

            self.Xbound = [-2*max(self.R_clad[0],self.R_clad[1]),2*max(self.R_clad[0],self.R_clad[1])]

            self.Ybound = [-max(self.R_clad[0],self.R_clad[1]),max(self.R_clad[0],self.R_clad[1])]

            self.cst = 2*(2-sqrt(2))

            self.delta = self.R_clad[0] + self.R_clad[1] - 2 * self.f * (self.R_clad[0] + self.R_clad[1] - np.sqrt(self.R_clad[0]**2 + self.R_clad[1]**2) )


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

        if self.OBJECT:
            [ ax.plot(*object.exterior.xy) for object in self.OBJECT]

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





        for obj in self.OBJ[:]:
            if obj.is_empty:
                continue
            patch = PolygonPatch(obj, alpha=0.5)
            ax1.add_patch(patch)



        ax2.set_title('Rasterized RI profil', fontsize=10)

        pcm = ax2.pcolormesh(np.linspace(*config.Y_bound, config.NY_point),
                             np.linspace(*config.X_bound, config.NX_point),
                             self.mesh,
                             cmap='PuBu_r',
                             norm=colors.LogNorm(vmin = 1.42, vmax=self.mesh.max()))


        if args.savefig:
            plt.savefig('result/Fusion{0:.2f}.png'.format(self.f))

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

        #self.OBJECT.append(self.overlap)

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

        self.OBJECT += [self.triangle_top, self.triangle_bottom]

        return Rv




    def construct_virtual(self, Rv):

        self.circ2, self.circ3  = self.C[0].buffer(Rv), self.C[1].buffer(Rv)

        self.OBJECT += [self.circ2, self.circ3]

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

        #self.OBJECT.append(self.added)


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

        self.OBJECT = [self.Clad[0], self.Clad[1]]

        Rv = self.construct_triangle_from_Rv(Rv)

        self.construct_virtual(Rv)

        self.create_mask()

        self.compute_added()

        cost = np.abs( self.added.area - self.overlap.area )

        if not args.quiet:
            sys.stdout.write('Topology:{0} \t Cost:{1:.5f} \t Rv:{2:.0f}\n'.format(self.topo, cost, float(Rv)))

        return cost




    def create_cores(self):

        self.OBJ = [self.coupler]

        N = 0

        for k, v in self.dict['geometry'].items():
            if k == 'Clad':
                continue

            if v['position'] == 'center':
                pos = (0,0)

            elif v['position'] == 'core0':
                pos = self.core_position[0]

            elif v['position'] == 'core1':
                pos = self.core_position[1]


            if v['function'] == 'clad':
                N += 1
                continue


            shape = Point(pos).buffer(v['radius'])

            self.OBJ.append(shape)

            N += 1



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

        if not args.quiet:
            sys.stdout.write('Cost:{0:.2f} \t position left:{1:.2f} \t position right:{2:.2f}\n'.format(cost0 + cost1, centers[0], centers[1]))

        return cost1, cost0


    def optimize_clad(self, init):


        if self.N_fiber == 1:
            self.coupler = Point([0,0]).buffer(self.R_clad[0])

            self.C = []
            self.P = []
            self.F = [Point([0,0])]

        if self.N_fiber == 2:

            self.construct_fibers()

            if self.topo == 'concave':
                bounds = ((self.R_clad[0]+self.delta+self.R_clad[1])/2,10000)
            else:
                bounds = (0,10000)

            z = least_squares(self.fusion_cost,init,bounds=bounds)




    def optimize_core_position(self, init=0):

        if self.N_fiber == 1:
            self.core_position = (0,0)

        if self.N_fiber == 2:
            z = fsolve(self.measure_side_area, [self.F[0].x, self.F[1].x])

            self.core_position = list( zip( z,[0,0] ) )



    def add_object_to_mesh(self, mesh, obj, index):

        mesh *= (-1 * obj + 1)

        mesh += obj * index

        return mesh



    def rasterize_polygone(self, polygone, points):

        obj_ext = Path(list( polygone.exterior.coords))

        obj_ext = obj_ext.contains_points(points).reshape((self.dict['mesh']['NX_point'],self.dict['mesh']['NY_point']))

        return obj_ext



    def rasterize_mesh(self):

        nx, ny = self.dict['mesh']['NX_point'], self.dict['mesh']['NY_point']

        self.mesh = np.ones((nx, ny))

        xv = np.linspace(*self.dict['mesh']['X_boundary'],nx)
        yv = np.linspace(*self.dict['mesh']['Y_boundary'],ny)

        x, y = np.meshgrid(xv,yv)

        x, y = x.flatten(), y.flatten()

        points = np.vstack((x,y)).T

        N = 0

        for k, v in self.dict['geometry'].items():

            if isinstance(v['index'], str):
                index = eval(v['index'])

            else:
                index = v['index']

            if args.debug:
                sys.stdout.write('Processing:{0} \t index:{1}'.format(k, index))

            Object = self.OBJ[N]

            obj = self.rasterize_polygone(Object, points)

            self.mesh = self.add_object_to_mesh(self.mesh, obj, index)

            N += 1


        if args.debug:
            self.plot_geometry()



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
