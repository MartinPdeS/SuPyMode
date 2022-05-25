#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import numpy as np
from descartes           import PolygonPatch
from shapely.geometry    import Point, LineString, MultiPolygon, Polygon

from SuPyMode.Plotting.PlotsUtils  import *
import matplotlib.pyplot as plt

try:
    from mayavi     import mlab
    from tvtk.tools import visual
except ImportError:
    logging.warning('Mayavi package could not be loaded! Not 3D rendering available.')


class Scene2D:
    UnitSize       = (4, 4)

    def __init__(self, nCols, nRows, ColorBar=True, Projection=None, Grid=True):
        self.nCols      = nCols
        self.nRows      = nRows
        self.ColorBar   = ColorBar
        self.Projection = Projection
        self.Grid       = Grid
        self.Boundaries = {'x': [0, 0], 'y': [0, 0]}

        self.InitScene()

    def InitScene(self):
        plt.rcParams['axes.grid'] = self.Grid
        plt.rcParams['ytick.labelsize'] = '8'
        plt.rcParams['xtick.labelsize'] = '8'
        plt.rcParams["font.size"]       = 10
        plt.rcParams["font.family"]     = "serif"

        FigSize = [ self.UnitSize[0]*self.nCols, self.UnitSize[1]*self.nRows ]

        self.Figure, self.Axes = plt.subplots(self.nRows,
                                              self.nCols,
                                              figsize=FigSize,
                                              subplot_kw=dict(projection=self.Projection))

        if np.size(self.Axes) == 1:
            self.Axes = np.reshape(np.asarray( [self.Axes] ), (1, 1))
        else:
            self.Axes = np.reshape(np.asarray( [self.Axes] ), (self.nRows, self.nCols))

    def AddLine(self, x, y, Col, Row, Title=None, Fill=True, Color='C0', xLabel=None, yLabel=None):
        ax = self.Axes[Row, Col]

        ax.plot(x, y, color='k')

        if Title:
            ax.set_title(Title)

        if Fill:
            ax.fill_between(x, y.min(), y, color=Color, alpha=0.7)

        if xLabel:
            ax.set_xlabel(xLabel)

        if yLabel:
            ax.set_ylabel(yLabel)


    def AddShapely(self, Col, Row, Object, Text=None):
        ax = self.Axes[Row, Col]

        if isinstance(Object, Point):
            ax.scatter(Object.x, Object.y, linewidth=7)
            self.UpdateBoundary(x=Object.x, y=Object.y)

        if isinstance(Object, (Polygon, MultiPolygon) ):
            Image = ax.add_patch( PolygonPatch(Object, alpha=0.5) )
            Coord = np.array( list(zip(*Object.exterior.coords.xy)) )
            self.UpdateBoundary(x=Coord[:,0], y=Coord[:,1])

        if Text:
            ax.text(Object.x, Object.y, Text)


    def UpdateBoundary(self, x=None, y=None):
        if x is not None:
            xMax, xMin = np.max(x), np.min(x)
            XBound = min(self.Boundaries['x'][0], xMin), max(self.Boundaries['x'][1], xMax)
            self.Boundaries['x'] = XBound

        if y is not None:
            yMax, yMin = np.max(y), np.min(y)
            YBound = min(self.Boundaries['y'][0], yMin), max(self.Boundaries['y'][1], yMax)
            self.Boundaries['y'] = XBound



    def AddMesh(self, Col, Row, x, y, Scalar, ColorMap='viridis', Title=None, xLabel=None, yLabel=None):
        ax = self.Axes[Row, Col]
        Image = ax.pcolormesh(x, y, Scalar, cmap=ColorMap, shading='auto')

        if Title is not None:
            ax.set_title(Title)

        if xLabel is not None:
            ax.set_xlabel(xLabel)

        if yLabel is not None:
            ax.set_ylabel(yLabel)

        if self.ColorBar:
            plt.colorbar(Image, ax=ax)

        self.UpdateBoundary(x=x, y=y)


    def SetLimits(self, Row, Col, XLim="Auto", YLim="auto"):
        ax = self.Axes[Row, Col]

        if XLim == 'auto': ax.set_xlim(self.Boundaries['x'])
        if YLim == 'auto': ax.set_ylim(self.Boundaries['y'])

        if isinstance(XLim, list): ax.set_xlim(XLim)
        if isinstance(YLim, list): ax.set_ylim(YLim)

    def Show(self):
        for ax in self.Axes.flatten():
            ax.set_aspect('equal')

        plt.tight_layout()
        plt.show()


class Scene3D:
    Size       = (600, 350)
    BackGround = (1,1,1)
    ForGround  = (0,0,0)

    def __init__(self, Axes=None, UnitSphere=None, ColorBar=True, AddSource=True, Bounds=None):
        self.Objects    = []
        self.Axes       = Axes
        self.UnitSphere = UnitSphere
        self.ColorBar   = ColorBar
        self.AddSource  = AddSource
        self.Bounds     = Bounds

        self.InitScene()


    def InitScene(self):
        self.Figure = mlab.figure(size    = self.Size,
                                  bgcolor = self.BackGround,
                                  fgcolor = self.ForGround
                                  )

        visual.set_viewer(self.Figure)


    def ConfigScene(self):
        self.Figure.scene.camera.position = [5, 8.5, 9.5]
        self.Figure.scene.camera.focal_point = [0.04, 0.8, -0.6]
        self.Figure.scene.camera.view_angle = 30.0
        self.Figure.scene.camera.view_up = [0.93, -0.22, -0.28]
        self.Figure.scene.camera.clipping_range = [0.03, 31]



    def AddStructured(self,
                      Coordinate,
                      Scalar,
                      Offsets         = [0,0,0],
                      ColorMap        = 'RdBu',
                      Opacity         = 1,
                      Name            = '',
                      Orientation     = 'horizontal',
                      Text            = None
                      ):

        Coordinate = [Coordinate[i] + Offsets[i] for i in range(3)]

        object = {'Type'           : 'Structured',
                  'Coordinate'     : Coordinate,
                  'Scalar'         : Scalar,
                  'Offsets'        : Offsets,
                  'ColorMap'       : ColorMap,
                  'Opacity'        : Opacity,
                  'Extra'          : True,
                  'Name'           : Name,
                  'Orientation'    : Orientation,
                  'Text'           : Text
                 }

        self.Objects.append( object )


    def AddQuiver(self,
                  Phi,
                  Theta,
                  Quiver,
                  Offsets = [0,0,0],
                  Skip    = 5,
                  Text    = None
                  ):

        Coordinate, PhiQuiver, ThetaQuiver = GetUnitVector(Phi[::Skip, ::Skip],
                                                           Theta[::Skip, ::Skip])

        Coordinate = [Coordinate[i]*1.1 + Offsets[i] for i in range(3)]

        if Quiver == 'Theta':
            Vector = ThetaQuiver
        elif Quiver == 'Phi':
            Vector = PhiQuiver


        object = {'Type'       : 'Quiver',
                  'Coordinate' : Coordinate,
                  'Vector'     : Vector,
                  'Offsets'    : Offsets,
                  'Color'      : (0,0,0),
                  'Scale'      : 0.1,
                  'Extra'      : False,
                  'Name'       : 'QuiverPlot',
                  'ScaleMode'  : 'vector',
                  'Rendered'   : False,
                  'Text'       : Text
                 }

        self.Objects.append( object )


    def AddUnstructured(self,
                        Coordinate,
                        Scalar,
                        Offsets         = [0,0,0],
                        ColorMap        = 'RdBu',
                        Opacity         = 1,
                        Name            = '',
                        Orientation     = 'horizontal',
                        Text            = None
                        ):

        Coordinate = [Coordinate[i] + Offsets[i] for i in range(3)]

        object = {'Type'       : 'Unstructured',
                  'Coordinate' : Coordinate,
                  'Scalar'     : Scalar,
                  'Offsets'    : Offsets,
                  'ColorMap'   : ColorMap,
                  'Opacity'    : Opacity,
                  'Extra'      : True,
                  'Name'       : Name,
                  'Orientation': Orientation,
                  'Rendered'   : False,
                  'Text'       : Text
                 }

        self.Objects.append( object )


    def AddExtra(self, Object):
        if self.Axes:
            AddUnitAxes(Figure=self.Figure, Scale=2.5, Origin=Object['Offsets'], ScaleTube=0.7)

        if self.UnitSphere:
            AddUnitSphere(Num = 50, Radius = 1., Origin = Object['Offsets'], Figure = self.Figure)

        if self.Bounds:
            Object['Image'].module_manager.scalar_lut_manager.data_range = self.Bounds

        if self.ColorBar and Object['ColorMap'] is not None:
            cb = mlab.colorbar(object            = Object['Image'],
                               label_fmt         = "%.0e",
                               nb_labels         = 5,
                               title             = Object['Name'],
                               orientation       = Object['Orientation'] )

            cb.scalar_bar.unconstrained_font_size = True
            cb.title_text_property.font_family    = 'times'
            cb.label_text_property.font_size      = 20

        if self.AddSource:
            Pad = [0, 0, -4]
            Offset = [Object['Offsets'][i] + Pad[i] for i in range(3)]
            WavenumberArrow(self.Figure, Origin=Offset, Scale=1)

        if Object['Text']:
            mlab.text3d(x          = Object['Offsets'][0],
                        y          = Object['Offsets'][1],
                        z          = Object['Offsets'][2]-2,
                        text       = Object['Text'],
                        line_width = 0.1,
                        figure     = self.Figure,
                        scale      = 0.25,
                        color      = (0,0,0))


    def RenderStructured(self, Object):

        Object['Image'] = mlab.mesh(*Object['Coordinate'],
                          scalars  = Object['Scalar'],
                          colormap = Object['ColorMap'],
                          figure   = self.Figure,
                          opacity  = Object['Opacity'] )

        Object['Rendered'] = True


    def RenderUnstructured(self, Object):
        Object['Image'] = mlab.points3d(*Object['Coordinate'], Object['Scalar'],
                                         mode       = 'sphere',
                                         scale_mode = 'none',
                                         colormap   = Object['ColorMap'])

        Object['Rendered'] = True


    def RenderQuiver(self, Object):
        Object['Image'] = mlab.quiver3d(*Object['Coordinate'],
                                        *Object['Vector'],
                                        color        = Object['Color'],
                                        scale_factor = Object['Scale'],
                                        scale_mode   = Object['ScaleMode'],
                                        figure       = self.Figure)

        Object['Rendered'] = True



    def Render(self):
        for Object in self.Objects:

            if Object['Type'] == 'Structured':
                self.RenderStructured(Object)


            elif Object['Type'] == 'Unstructured':
                self.RenderUnstructured(Object)

            elif Object['Type'] == 'Quiver':
                self.RenderQuiver(Object)

            if Object['Extra']:
                self.AddExtra(Object)


    def Show(self, Wait=False):
        self.ConfigScene()

        if Wait: return

        mlab.show()




# -
