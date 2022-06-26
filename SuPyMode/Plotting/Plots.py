#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import numpy as np
from descartes           import PolygonPatch
from shapely.geometry    import Point, LineString, MultiPolygon, Polygon
import matplotlib.pyplot as plt
import matplotlib.cm     as cm
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
    def __init__(self, X, Y, Scalar, ColorMap='viridis', Label=''):
        self.X = X
        self.Y = Y
        self.Scalar = Scalar
        self.ColorMap=ColorMap
        self.Label = Label


    def Render(self, Ax, ColorBarDict):
        self.Norm = colors.BoundaryNorm(np.unique(self.Scalar), 200, extend='both') if ColorBarDict['Discreet'] else None

        Image = Ax._ax.pcolormesh(self.X,
                              self.Y,
                              self.Scalar,
                              cmap    = self.ColorMap,
                              shading = 'auto',
                              norm = self.Norm
                              )

        if ColorBarDict['Create']:
            divider = make_axes_locatable(Ax._ax)
            cax = divider.append_axes(ColorBarDict['Position'], size="10%", pad=0.15)

            if ColorBarDict['Discreet']:
                ticks = np.unique(self.Scalar)
                plt.colorbar(mappable=Image, norm=self.Norm, boundaries=ticks, ticks=ticks, cax=cax)
            else:
                plt.colorbar(mappable=Image, norm=self.Norm, cax=cax)


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
    def __init__(self,
                 Row, Col,
                 xLabel='', yLabel='', Title='',
                 Grid=True, Legend=True,
                 xScale='linear', yScale='linear',
                 xLimits=None, yLimits=None, Equal=False,
                 ColorBar=False, ColorbarPosition='bottom', DiscreetColorbar=False, ):
        self.Row     = Row
        self.Col     = Col
        self.Legend  = Legend
        self.Artist  = []
        self._ax     = None
        self.Grid    = Grid
        self.xScale  = xScale
        self.yScale  = yScale
        self.xLimits = xLimits
        self.yLimits = yLimits
        self.Equal   = Equal

        self.Labels  = {'x': xLabel,
                        'y': yLabel,
                        'Title': Title}

        self.ColorBarDict = {'Create': ColorBar,
                             'Position': ColorbarPosition,
                             'Discreet': DiscreetColorbar}


    def AddArtist(self, *Artist):
        for art in Artist:
            self.Artist.append(art)

    def Render(self):
        for art in self.Artist:
            Image = art.Render(self, ColorBarDict=self.ColorBarDict)

        if self.Legend:
            self._ax.legend()

        self._ax.grid(self.Grid)

        if self.xLimits is not None: self._ax.set_xlim(self.xLimits)
        if self.yLimits is not None: self._ax.set_ylim(self.yLimits)

        self._ax.set_xlabel(self.Labels['x'])
        self._ax.set_ylabel(self.Labels['y'])
        self._ax.set_title(self.Labels['Title'])

        self._ax.set_xscale(self.xScale)
        self._ax.set_yscale(self.yScale)

        if self.Equal:
            self._ax.set_aspect("equal")






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
        if Ax.ndim == 1: Ax = np.asarray([Ax])

        self.Figure.suptitle(self.Title)

        for ax in self.Axis:
            ax._ax = Ax[ax.Row, ax.Col]


    def Render(self):
        self.GenerateAxis()

        for ax in self.Axis:
            ax.Render()

        plt.tight_layout()


    def Show(self):
        self.Render()
        plt.show()






# -
