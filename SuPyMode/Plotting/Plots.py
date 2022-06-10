#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import numpy as np
from descartes           import PolygonPatch
from shapely.geometry    import Point, LineString, MultiPolygon, Polygon
import matplotlib.pyplot as plt
from matplotlib          import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

from SuPyMode.Plotting.PlotsUtils  import FieldMap, MidPointNorm
from SuPyMode.Tools.utils import ToList


try:
    from mayavi     import mlab
    from tvtk.tools import visual
except ImportError:
    logging.warning('Mayavi package could not be loaded! Not 3D rendering available.')

import matplotlib
matplotlib.style.use('ggplot')





class AddShapely:
    def __init__(self, Object, Text=None):
        self.Object = Object
        self.Text   = Text

    def Render(self, Ax):
        if isinstance(Object, Point):
            Ax.scatter(Object.x, Object.y, linewidth=7)

        if isinstance(Object, (Polygon, MultiPolygon) ):
            Image = Ax.add_patch( PolygonPatch(Object, alpha=0.5) )

        if Text:
            Ax.text(Object.x, Object.y, Text)





class Contour:
    def __init__(self, X, Y, Scalar, ColorMap='viridis', Title=None, xLabel=None, yLabel=None, IsoLines=None):
        self.X = X
        self.Y = Y
        self.Scalar = Scalar
        self.ColorMap = ColorMap
        self.Label = Label
        self.IsoLines = IsoLines


    def Render(self, Ax):
        Image = Ax.contour(self.X,
                            self.Y,
                            self.Scalar,
                            level = self.IsoLines,
                            colors="black",
                            linewidth=.5 )

        Image = Ax.contourf(self.X,
                            self.Y,
                            self.Scalar,
                            level = self.IsoLines,
                            cmap=self.ColorMap,
                            norm=colors.LogNorm() )



class Mesh:
    def __init__(self, X, Y, Scalar, ColorMap='viridis', DiscretNorm=None, Label=''):
        self.X = X
        self.Y = Y
        self.Scalar = Scalar
        self.ColorMap=ColorMap
        self.Label = Label

        #norm = colors.BoundaryNorm(DiscretNorm, ColorMap) if DiscretNorm is not None else None

    def Render(self, Ax):
        Image = Ax.pcolormesh(self.X,
                              self.Y,
                              self.Scalar,
                              cmap=self.ColorMap,
                              shading='auto',
                              vmin=-np.max(np.abs(self.Scalar)),
                              vmax=+np.max(np.abs(self.Scalar)))

        return Image



class Line:
    def __init__(self, X, Y, Label=None, Fill=False, Color=None):
        self.X = X
        self.Y = Y
        self.Fill = Fill
        self.Label = Label
        self.Color = Color

    def Render(self, Ax):

        Ax.plot(self.X, self.Y, label=self.Label)

        if self.Fill:
            Ax.fill_between(self.X, self.Y.min(), self.Y, color=self.Color, alpha=0.7)



class Axis:
    def __init__(self, Row, Col, xLabel, yLabel, Title, Grid=True, Legend=True, xScale='linear', yScale='linear', xLimits=None, yLimits=None, Equal=False, ColorBar=False, ColorbarPosition='bottom'):
        self.Row    = Row
        self.Col    = Col
        self.xLabel = xLabel
        self.yLabel = yLabel
        self.Title  = Title
        self.Legend = Legend
        self.Artist = []
        self.MPLAxis = None
        self.Grid    = Grid
        self.xScale  = xScale
        self.yScale  = yScale
        self.xLimits = xLimits
        self.yLimits = yLimits
        self.Equal   = Equal
        self.ColorBar = ColorBar
        self.ColorbarPosition = ColorbarPosition


    def AddArtist(self, *Artist):
        for art in Artist:
            self.Artist.append(art)

    def Render(self):
        for art in self.Artist:
            Image = art.Render(self.MPLAxis)

        if self.Legend:
            self.MPLAxis.legend()

        self.MPLAxis.grid(self.Grid)

        if self.xLimits is not None: self.MPLAxis.set_xlim(self.xLimits)
        if self.yLimits is not None: self.MPLAxis.set_ylim(self.yLimits)

        self.MPLAxis.set_xlabel(self.xLabel)
        self.MPLAxis.set_ylabel(self.yLabel)
        self.MPLAxis.set_title(self.Title)

        self.MPLAxis.set_xscale(self.xScale)
        self.MPLAxis.set_yscale(self.yScale)

        if self.Equal:
            self.MPLAxis.set_aspect("equal")

        if self.ColorBar:
            divider = make_axes_locatable(self.MPLAxis)
            cax = divider.append_axes(self.ColorbarPosition, size="5%", pad=0.05)
            plt.colorbar(Image, cax=cax, )




class Scene:
    UnitSize = (10, 3)
    plt.rcParams['ytick.labelsize'] = 8
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams["font.size"]       = 8
    plt.rcParams["font.family"]     = "serif"
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['axes.linewidth'] = 1.5

    def __init__(self, Title='', UnitSize=None):
        self.Axis = []
        self.Title = Title
        self.nCols = 1
        self.nRows = None
        if UnitSize is not None: self.UnitSize = UnitSize


    def AddAxes(self, *Axis):
        for ax in Axis:
            self.Axis.append(ax)


    def GetMaxColsRows(self):
        RowMax, ColMax = 0,0
        for ax in self.Axis:
            RowMax = ax.Row if ax.Row > RowMax else RowMax
            ColMax = ax.Col if ax.Col > ColMax else ColMax

        return RowMax, ColMax


    def GenerateAxis(self):
        RowMax, ColMax = self.GetMaxColsRows()

        self.nRows = len(self.Axis)

        FigSize = [ self.UnitSize[0]*(ColMax+1), self.UnitSize[1]*(RowMax+1) ]

        self.Figure, Ax  = plt.subplots(ncols=ColMax+1, nrows=RowMax+1, figsize=FigSize)

        if not isinstance(Ax, np.ndarray): Ax = np.asarray([[Ax]])
        if Ax.shape == 1: Ax = np.asarray([A])

        self.Figure.suptitle(self.Title)

        for ax in self.Axis:
            ax.MPLAxis = Ax[ax.Row, ax.Col]


    def Render(self):
        self.GenerateAxis()

        for ax in self.Axis:
            ax.Render()

        plt.tight_layout()


    def Show(self):
        self.Render()
        plt.show()

















class Scene2D:
    UnitSize       = (10, 3)

    def __init__(self, nCols, nRows, ColorBar=True, Projection=None, Grid=True, UnitSize=None):
        self.nCols      = nCols
        self.nRows      = nRows
        self.ColorBar   = ColorBar
        self.Projection = Projection
        self.Grid       = Grid
        self.Boundaries = {'x': [0, 0], 'y': [0, 0]}
        if UnitSize is not None: self.UnitSize = UnitSize

        self.InitScene()

    def InitScene(self):
        plt.rcParams['axes.grid'] = self.Grid
        plt.rcParams['ytick.labelsize'] = 8
        plt.rcParams['xtick.labelsize'] = 8
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


    def AddLine(self, x, y, Col, Row, Title=None, Fill=True, Color=None, xLabel=None, yLabel=None, Legend=None):
        ax = self.Axes[Row, Col]

        ax.plot(x, y, color=Color, label=Legend)

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


    def SetAxes(self, Col, Row, Equal=None, Legend=None, xLimits=None, yLimits=None, xScale='linear', yScale='linear', LegendTitle=None):
        ax = self.Axes[Row, Col]

        if Equal is True:
            ax.set_aspect('equal')
        if Equal is False:
            ax.set_aspect('auto')

        if Legend is not None:
            ax.legend().set_visible(Legend)

        if xLimits == 'auto': ax.set_xlim(self.Boundaries['x'])
        if yLimits == 'auto': ax.set_ylim(self.Boundaries['y'])

        if isinstance(xLimits, list): ax.set_xlim(xLimits)
        if isinstance(yLimits, list): ax.set_ylim(yLimits)

        ax.set_yscale(yScale)
        ax.set_xscale(xScale)

        if LegendTitle: ax.legend().set_title(LegendTitle)



    def AddMesh(self, Col, Row, x, y, Scalar, ColorMap='viridis', Title=None, xLabel=None, yLabel=None, DiscretNorm=None):
        ax = self.Axes[Row, Col]

        norm = colors.BoundaryNorm(DiscretNorm, ColorMap) if DiscretNorm is not None else None

        Image = ax.pcolormesh(x, y, Scalar, cmap=ColorMap, shading='auto', vmin=-np.max(np.abs(Scalar)), vmax=+np.max(np.abs(Scalar)) )

        if Title is not None:
            ax.set_title(Title)

        if xLabel is not None:
            ax.set_xlabel(xLabel)

        if yLabel is not None:
            ax.set_ylabel(yLabel)

        if self.ColorBar:
            plt.colorbar(Image, ax=ax, location='bottom')

        self.UpdateBoundary(x=x, y=y)
        plt.tight_layout(pad=3)


    def AddContour(self, Col, Row, x, y, Scalar, ColorMap='viridis', Title=None, xLabel=None, yLabel=None, IsoLines=None):
        ax = self.Axes[Row, Col]

        ax.contour(x, y, Scalar, levels=IsoLines, colors='black', linewidth=.5)
        Image = ax.contourf(x, y, Scalar, levels=IsoLines, cmap='jet', norm=colors.LogNorm() )

        if Title is not None:
            ax.set_title(Title)

        if xLabel is not None:
            ax.set_xlabel(xLabel)

        if yLabel is not None:
            ax.set_ylabel(yLabel)

        if self.ColorBar:
            plt.colorbar(Image, ax=ax, location='bottom', format="%.3f")

        self.UpdateBoundary(x=x, y=y)
        plt.tight_layout(pad=3)


    def Show(self):
        #plt.tight_layout(pad=3)
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
