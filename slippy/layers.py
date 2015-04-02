#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Layers class of slippy.

:copyright:
    Dan Sandiford (sonderfjord@gmail.com), 2015
    Louis Moresi, 2015
:license:
    MIT
"""
import math
import numpy as np
from shapely.affinity import translate
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.geometry import Point
import underworld as uw

def create_layers(layer_depths,  trench_shp, Rc = 180, angle = 70, Zmax = 250, lengthscale = 111000., zmax = 10, zmin = 0, slippy_width = 20):
    """
    Functionto that creates functions representing slabs at depth, 
    based on radius of curvature, slope at which the slab becomes straight
    maximum depth (of the top surface), etc. 
    Layer depths are given in kilometers. It returns depths in model units.
    The x, y coordinates are assumed to be in degrees (lat, lon), but are used in a cartesian sense.
    Hence, it's best if the model is situated near (0,0)

    
    This function will create layers that extend from the trench in the "y" direction. 
    It is assumed that the reference frame has been shifted so that relative plate motion is North-South, 
    i.e. the Y axis in python/numpy 
    
    Depths are considered positive
    """
    from scipy import interpolate
    #Convert to model coordinates:
    mRc = Rc*1000./lengthscale
    def create_surface(ldepth):
        #Within this additional function model-space coordinates are assumed
        rl = mRc - ldepth  #adjusted radius of curvature.
        ang_rad = angle*(math.pi/180)
        slope = math.tan(ang_rad)
        if slope > 0:
            slope = slope*-1.
        else:
            slope = slope
        z_ec = math.sqrt((rl**2/ ((slope**2)+1)))
        y_ec = math.sqrt(rl**2 - z_ec**2)            
        """
        Here we define the depth function f(x,y)
        """
        def surface(x,y):
            tpts = list(trench_shp.coords)
            px = [i[0] for i in tpts]
            py = [i[1] for i in tpts]
            f = interpolate.interp1d(px, py, kind='linear', axis=-1)
            px = x
            YaC = abs(y - float(f(px)))            #Model units
            #print "YaC is %s." % YaC
            if YaC < y_ec:              
                #Ykm = YaC*111.
                ZaC = math.sqrt(rl**2 - YaC**2) #absolute height from refernce depth of rl
                if ZaC > mRc:                   #Function is zero if above surface
                    fdepth = mRc
                else:
                    fdepth = ZaC                #absolute height from refernce depth of rl

            else:
                YaS = abs(y - float(f(px))) - y_ec    #distance from the curved - straight crossover point
                fdepth = z_ec + (YaS * slope)     #not sure about ldepth          
            return fdepth + zmax - mRc
        return surface
        
    """
    Make the layers:
    """
    fdepths = list(layer_depths) #Deep copy
    #Add the zero depth layer
    fdepths.insert(0, 0.)
    #Add the Slippy layer
    fdepths.insert(0, float(-1*slippy_width))
    if len(fdepths) == 2:
        layer_depths.reverse
    mdep = []
    [mdep.append(i*1000./lengthscale) for i in fdepths]
    f_basket = []
    for dep in mdep:
        f = create_surface(dep)
        f_basket.append(f)
    return f_basket
    
       
    
def create_slab_region(trench_shp,direction = 1, Rc = 180, angle = 70, Zmax = 250, lengthscale = 111000.):
    """
    This function creates creates a polygon that covers the region
    where the slab will be interpolated, i.e. from the trench to the point
    defined by the parameters: 
    Radius of curvature, 
    Angle after which the slab is flat
    Zmax: the deepest part of the (top) slab
    """
    mZ = Zmax*1000./lengthscale
    mRc = Rc*1000./lengthscale
    ang_rad = angle*(math.pi/180)
    slope = math.tan(ang_rad)
    if slope > 0:
        slope = slope*-1.
    else:
        slope = slope
    #print slope
    z_ec = math.sqrt((mRc**2/ ((slope**2)+1)))
    y_ec = math.sqrt(mRc**2 - z_ec**2)
    
    #work out northward distance covering the straight part of the slab
    if mZ >= z_ec:
        tdep = mZ - (mRc - z_ec)
    else:
        tdep = 0
    y2 = abs(tdep / slope)
    tdis = y2 + y_ec                      #model units
    if direction == 1:
        pass
    elif direction == -1:
        tdis = tdis*-1.
    print tdis
    #Move trench northward by y. 
    trench = trench_shp
    #trench1 = translate(trench, xoff=0.0, yoff=my, zoff=0.0)
    trench2 = translate(trench, xoff=0.0, yoff=tdis, zoff=0.0)
    tpts = list(trench.coords)
    #Create shapefiles for the curved and straight slab part
    tpts2 = list(trench2.coords)
    tpts2.reverse()
    temp_crds = tpts + tpts2
    slab_region = Polygon(temp_crds)
    return(slab_region)
    
    
    
    
def map_flat_layers(position, dim, lon = [], lat = [], shapes=[], layers = [], upper_mantle_index=0, lower_mantle_index=1, lower_mantle_depth = 4, mantle_buffer = [10,-10,10,-10]):
    """
    This function loops through shapes (Shapely Polygons) and layers, and assigns material indexes.
    """
    if dim == 3:
        lat = position[1]
        lon = position[0]
        depth = position[2]
    elif dim == 2:
        if lon:
            lon = lon[0]
            lat = position[0]
            depth = position[1]
        elif lat:
            lat = lat[0]
            lon = position[0]
            depth = position[1]
        else:
            cs = raw_input("enter x or y axis to make cross-section: ")
            if cs == "x":
                lon = position[0]
                lat = raw_input("enter y-pos to make x-cross-section: ")
                lat = float(lon)
                depth = position[1]
            elif cs == "y":
                lat = position[0]
                lon = raw_input("enter x-pos to take y-cross-section")
                lon = float(lat)
                depth = position[1]
            else:
                print "Must enter x or y"
    else:
        print "dim must be 2 or 3"
    #any points outside mantle buffer become upper mantle   
    #any point beneath upper mantle become lower mantle
    if depth < lower_mantle_depth:
        mat_index = lower_mantle_index     
    elif lon > mantle_buffer[0]:
        mat_index = upper_mantle_index
    elif lon < mantle_buffer[1]:
        mat_index = upper_mantle_index
    elif lat > mantle_buffer[2]:
        mat_index = upper_mantle_index
    elif lat < mantle_buffer[3]:
        mat_index = upper_mantle_index
    else:
    #use shapely point class
    #pyGplates to be used in near-future
        pt = Point(lon, lat)
        found_shape = 0
        mat_index = 0
        mat_counter = 0
        for i in range(0, len(shapes)):
            if i == 0:
                mat_counter = 4
            else: 
                mat_counter = mat_counter + ((len(layers[i-1])-1))
            mat_index = mat_counter
            found_layer = 0
            if found_shape ==1:
                break
            elif shapes[i].contains(pt):
                found_shape == 1
                found_layer = 0 #reset after layer loop
                for j in range(0, len(layers[i][:-1])):
                    k = j+1
                    if found_layer == 1:
                        break
                    elif layers[i][j] > depth > layers[i][k]:
                        mat_index = mat_index
                        found_layer = 1
                        break
                    elif k == len(layers[i]) -1:
                        mat_index = upper_mantle_index
                        break
                    else:
                        mat_index = mat_index +1    
                break    
            else:
                mat_index = upper_mantle_index
    return mat_index   
    
    
    
def map_slab_layers(ii, mat_indx_var, position, dim, cutoff  = 250, lon = [], lat = [], shapes=[], layers = [], mantle_buffer = [10,-10,10,-10], stable_indx = 2., slippy_indx = 3.):
    """
    This function loops through shapes (Shapely Polygons) and layers, and assigns material indexes.
    """
    if dim == 3:
        lat = position[1]
        lon = position[0]
        depth = position[2]
    elif dim == 2:
        if lon:
            lon = lon[0]
            lat = position[0]
            depth = position[1]
        elif lat:
            lat = lat[0]
            lon = position[0]
            depth = position[1]
        else:
            cs = raw_input("enter x or y axis to make cross-section: ")
            if cs == "x":
                lon = position[0]
                lat = raw_input("enter y-pos to make x-cross-section: ")
                lat = float(lon)
                depth = position[1]
            elif cs == "y":
                lat = position[0]
                lon = raw_input("enter x-pos to take y-cross-section")
                lon = float(lat)
                depth = position[1]
            else:
                print "Must enter x or y"
    else:
        print "dim must be 2 or 3"
    #No slab above Z max    
    if depth < cutoff:
        mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0])    
    elif lon > mantle_buffer[0]:
        mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0]) 
    elif lon < mantle_buffer[1]:
        mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0]) 
    elif lat > mantle_buffer[2]:
        mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0]) 
    elif lat < mantle_buffer[3]:
        mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0])
    #use shapely point class
    else: 
        pt = Point(lon, lat)
        found_shape = 0
        mat_index = 0
        mat_counter = 0
        for i in range(0, len(shapes)):
            #Mat counter starts at zero and adds the number of layers corresponding to each shape loop...
            mat_counter = mat_counter + ((len(layers[i])-2)) 
            mat_index = mat_counter
            found_layer = 0       
            if found_shape ==1:
                break
            elif shapes[i].contains(pt):
                found_shape == 1
                #First check if the slab layers are intended as "stabilizing regions"
                if len(layers[i])==2:
                    if layers[i][0](lon,lat) >= depth > layers[i][1](lon,lat):
                        mat_index = stable_indx
                        break
                    else:
                        mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0])   
                        break               
                else:
                    found_layer = 0 #reset after layer loop
                    for j in range(0, len(layers[i][:-1])):   #This -2 thing would change depending on how many layers 
                        k = j+1
                        if found_layer == 1:
                            break
                        elif layers[i][j](lon,lat) >= depth > layers[i][k](lon,lat):
                            #print(layers[i][j](lon,lat))
                            #The first layer is always the slippy layer, assign appropriate index, dont change mat_index
                            if j == 0:
                                mat_index = slippy_indx
                                found_layer = 1
                                break
                            else:
                                mat_index = mat_index
                                found_layer = 1
                                break
                        elif k == len(layers[i]) -1:
                            mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0])
                            found_layer = 1
                            break
                        else:
                            #Only increase mat_index once past "slippy" layer
                            if j == 0:            
                                mat_index = mat_index
                            else:
                                mat_index = mat_index +1   
                    break 
            #This is basically "do nothing return go to next item in shapes"
            else:
                mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0])
            
    return int(mat_index)