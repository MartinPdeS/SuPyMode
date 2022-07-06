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
from dataclasses import dataclass




try:
    from mayavi     import mlab
    from tvtk.tools import visual
except ImportError:
    logging.warning('Mayavi package could not be loaded! Not 3D rendering available.')

import matplotlib
matplotlib.style.use('ggplot')


@dataclass
class ColorBar:
    Color: str = 'viridis'
    Discreet: bool = False
    Position: str = 'left'
    Orientation   = "vertical"


@dataclass
class AddShapely:
    Object: list
    Text: str = None

    def Render(self, Ax):
        if isinstance(Object, Point):
            Ax.scatter(Object.x, Object.y, linewidth=7)

        if isinstance(Object, (Polygon, MultiPolygon) ):
            Image = Ax.add_patch( PolygonPatch(Object, alpha=0.5) )

        if Text:
            Ax.text(Object.x, Object.y, Text)




@dataclass
class Contour:
    X: np.ndarray
    Y: np.ndarray
    Scalar: np.ndarray
    ColorMap: str = 'viridis'
    xLabel: str = ''
    yLabel: str = ''
    IsoLines: list = None

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


@dataclass
class Mesh:
    X: np.ndarray
    Y: np.ndarray
    Scalar: np.ndarray
    ColorMap: str = 'viridis'
    Label: str = ''

    def Render(self, Ax):
        self.Norm = colors.BoundaryNorm(np.unique(self.Scalar), 200, extend='both') if Ax.Colorbar.Discreet else None

        Image = Ax._ax.pcolormesh(self.X,
                              self.Y,
                              self.Scalar,
                              cmap    = self.ColorMap,
                              shading = 'auto',
                              norm = self.Norm
                              )

        if Ax.Colorbar is not None:
            divider = make_axes_locatable(Ax._ax)
            cax = divider.append_axes(Ax.Colorbar.Position, size="10%", pad=0.15)

            if Ax.Colorbar.Discreet:
                ticks = np.unique(self.Scalar)
                plt.colorbar(mappable=Image, norm=self.Norm, boundaries=ticks, ticks=ticks, cax=cax, orientation=Ax.Colorbar.Orientation)
            else:
                plt.colorbar(mappable=Image, norm=self.Norm, cax=cax, orientation=Ax.Colorbar.Orientation)


        return Image


@dataclass
class Line:
    X: np.ndarray
    Y: np.ndarray
    Label: str = None
    Fill: bool = False
    Color: str = None

    def Render(self, Ax):

        Ax._ax.plot(self.X, self.Y, label=self.Label)

        if self.Fill:
            Ax._ax.fill_between(self.X, self.Y.min(), self.Y, color=self.Color, alpha=0.7)



@dataclass
class Axis:
    Row: int
    Col: int
    xLabel: str = ''
    yLabel: str = ''
    Title: str = ''
    Grid: bool = True
    Legend: bool = True
    xScale: str = 'linear'
    yScale: str = 'linear'
    xLimits: list = None
    yLimits: list = None
    Equal: bool = False
    Colorbar: ColorBar = None

    def __post_init__(self):

        self._ax = None
        self.Artist  = []

        self.Labels  = {'x': self.xLabel,
                        'y': self.yLabel,
                        'Title': self.Title}


    def AddArtist(self, *Artist):
        for art in Artist:
            self.Artist.append(art)

    def Render(self):
        for art in self.Artist:
            Image = art.Render(self)

        if self.Legend:
            self._ax.legend(fancybox=True, facecolor='white', edgecolor='k')

        self._ax.grid(self.Grid)

        if self.xLimits is not None: self._ax.set_xlim(self.xLimits)
        if self.yLimits is not None: self._ax.set_ylim(self.yLimits)

        self._ax.set_xlabel(self.Labels['x'])
        self._ax.set_ylabel(self.Labels['y'])
        self._ax.set_title(self.Labels['Title'])

        self._ax.set_xscale(self.xScale)
        self._ax.set_yscale(self.yScale)

        self._ax.text(0.5, 0.1, 'SuPyMode', transform=self._ax.transAxes,
                fontsize=30, color='white', alpha=0.2,
                ha='center', va='baseline', rotation='0')

        if self.Equal:
            self._ax.set_aspect("equal")






class Scene:
    UnitSize = (10, 3)
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams["font.size"]       = 12
    plt.rcParams["font.family"]     = "serif"
    plt.rcParams['axes.edgecolor']  = 'black'
    plt.rcParams['axes.linewidth']  = 1.5

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
        print(Ax.shape)
        for ax in self.Axis:
            ax._ax = Ax[ax.Col, ax.Row]


    def Render(self):
        self.GenerateAxis()

        for ax in self.Axis:
            ax.Render()

        plt.tight_layout()


    def Show(self):
        self.Render()
        plt.show()






# -
