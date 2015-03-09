#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Some utility functionality.

:copyright:
    Dan Sandiford (sonderfjord@gmail.com), 2015
    Louis Moresi, 2015
:license:
    MIT
"""

import math
import numpy as np

###################
#Rheology and temperture relationships
###################

# rheology & temperature profile functions
def viscosity(T,D, E0, V0, R0):
    return math.exp((E0 + V0*D)/(R0*T))
    

def hsct(depth,time, Tint, Tsurf, kappa):
    Th = Tsurf +(Tint-Tsurf) * math.erf(depth/(2*math.sqrt(time*kappa)))
    return Th
    
def linear_geotherm(depth, Zthick, Tint, Tsurf):
    Tl = (Tint-Tsurf)/float(Zthick)
    Tl = Tl*depth
    return Tl
    
###################
#General utils
###################
    
    
def layers_rescale(length_scale,layers, max_model_vert=10):
    """
    Rescale layers or slab function to uw coordinates
    """
    scaled_layers = []
    if type(layers[0]) == int:
        for i in layers:
            depth = 1000.*i/length_scale
            out = -1*depth + max_model_vert
            scaled_layers.append(out)
    else:
        pass
    return scaled_layers 