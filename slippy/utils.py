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
import underworld as uw
import underworld.shapes as shape

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
    
def swarm_fom_box(uwdict, Shapely, num, name="dummy_name"):
    """
    This function takes a shapely box and a swarm,
    used for 2d swarms on the surface
    Needs to be extended to arbitary polygons
    
    """
    box_name = name + "_box"
    shape.setup.boxCreate(componentName=box_name , startList=[Shapely.bounds[0],Shapely.bounds[1],9.5], endList=[Shapely.bounds[2],Shapely.bounds[3],10])
    layout_name = name + "_layout"
    uwdict["components"][layout_name] = dict([('Type', 'WithinShapeParticleLayout'), ('shape', box_name ), ('totalInitialParticles', num)])
    swarm_name = name + "_swarm"
    uw.swarms.tools.TracerSwarmCreate(particleLayout=None, componentName=swarm_name, meshName='linearMesh')
    uwdict["components"][swarm_name]['ParticleLayout'] = layout_name
    uw.dictionary.SetDictionary(uwdict)