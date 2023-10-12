#!/usr/bin/env python
# -*- coding: utf-8 -*-

beta = {
    "show_legend": True,
    "x_label": 'Inverse taper ratio',
    "y_label": 'Propagation constant [rad/M]',
    "y_scale": "linear",
    "line_width": 2
}

index = {
    "show_legend": True,
    "x_label": 'Inverse taper ratio',
    "y_label": 'Effective refraction index',
    "y_scale": "linear",
    "y_limits": [1.44, 1.455],
    "line_width": 2
}

eigen_value = {
    "show_legend": True,
    "x_label": 'Inverse taper ratio',
    "y_label": 'Mode eigen values',
    "y_scale": "linear",
    "line_width": 2
}

field = {
    "show_legend": False,
    "x_label": r'X-Direction [$\mu m$]',
    "y_label": r'Y-direction [$\mu m$]',
    'x_scale_factor': 1e6,
    'y_scale_factor': 1e6,
    "equal": True
}

normalized_coupling = {
    "show_legend": True,
    "x_label": 'Inverse taper ratio',
    "y_label": 'Mode coupling',
    "y_scale": "linear",
    "line_width": 2
}

overlap = {
    "show_legend": True,
    "x_label": 'Inverse taper ratio',
    "y_label": 'Mode overlap intergral',
    "y_scale": "linear",
    "line_width": 2
}

beating_length = {
    "show_legend": True,
    "x_label": 'Inverse taper ratio',
    "y_label": 'Beating length [m]',
    "y_scale": "log",
    "line_width": 2
}

adiabatic = {
    "show_legend": False,
    "x_label": 'Inverse taper ratio',
    "y_label": r'Adiabatic criterion [$\mu$m$^{-1}$]',
    "y_scale": 'log',
    "y_scale_factor": 1e-6,
    "y_limits": [1e-5, 1],
    "line_width": 2
}

z_profile = {
    "show_legend": False,
    "x_label": 'Z-propagation [mm]',
    "y_label": 'Inverse taper ratio',
    "x_scale_factor": 1e3,
    "y_scale": "linear",
    "line_width": 2
}

taper_angle = {
    "show_legend": False,
    "y_label": 'Taper angle [rad]',
    "x_label": 'Z-propagation [mm]',
    "x_scale_factor": 1e3,
    "y_scale": "linear",
    "line_width": 2
}


# -
