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
from slippy import rotations

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
    
def midswarm(minX, maxX, minY, maxY, shapes=[], num=10000, depth = 9.545):
    """
    This defines a swarm at a given depth across the entire domain, 
    then separates into individual domains. 
    The way the model domain is passed is in very ugly, this needs to be standardised in V.2
    """
    import numpy as np
    from shapely.geometry import Point
    #####
    #Create the initial swarm
    #####
    xs = np.random.uniform(minX,maxX, num)
    ys = np.random.uniform(minY,maxY, num)
    zs = np.array([num, 1])
    zs = np.zeros(num)
    zs.fill(depth)
    mswarm = np.column_stack((xs,ys,zs))
    ####
    #Separate swarm on shapes
    ####
    final_swarms = []
    for pol in shapes:
        indx_list = []
        for i in range(0, len(mswarm[:,0])):
            pt = Point((mswarm[i,0],mswarm[i,1]))
            if pol.contains(pt):
                indx_list.append(i)
        swarm = mswarm[indx_list,:]
        final_swarms.append(swarm)
    return final_swarms
    
    
    
def rotate_grid_points(lon, lat, data, rm):
    '''
    A function to rotate geographic data, 
    given a rotation matrix rm. Works for arbitarry points, 
    as well as gridded data.
    '''
    #rotations.rotate_lat_lon(lat[n], lon[n], ra, ang)
    #convert rm to euler:
    ra, ang = rotations._get_axis_and_angle_from_rotation_matrix(rm)
    #first the case of arbitrary points
    if len(lon.shape) == len(lat.shape) == 1:
        nlons = []
        nlats = []
        for n in range(0, len(lat)):
            temp = rotations.rotate_lat_lon(lat[n], lon[n], ra, ang)
            nlons.append(temp[1])
            nlats.append(temp[0])
        alons = np.asarray(nlons)
        alats = np.asarray(nlats)        
        out = np.column_stack((alons, alats, data))
        return out
    elif len(lon.shape) == len(lat.shape) > 1:
        nlons = []
        nlats = []
        for i in range(0, lat.shape[0]):
            for j in range(0, lon.shape[1]):
                temp = rotations.rotate_lat_lon(lat[i,j], lon[i,j], ra, ang)
                nlons.append(temp[1])
                nlats.append(temp[0])
        alons = np.asarray(nlons)
        alats = np.asarray(nlats)
        flons = alons.reshape(lon.shape)
        flats = alats.reshape(lon.shape)
        out = np.dstack((flons, flats, data))
        return out
        
def spatial_mask(minX, maxX, minY, maxY, xy_array):
    """
    Simple function to get rid of values outside bounding box.
    No error checking...
    """
    indx_list = []
    for i in range(0,len(xy_array[:,0])):
        if minX < xy_array[i,0] < maxX and minY < xy_array[i,1] < maxY:
            indx_list.append(i)
    out = xy_array[indx_list,:]
    return out  
        