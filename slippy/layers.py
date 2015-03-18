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

def create_layers(layer_depths,  trench_shp, Rc = 180, angle = 70, Zmax = 250, kms_to_model = 111. ):
    """
    This function creates functions representing slabs at depth, 
    based on radius of curvature, slope at which the slab becomes straight
    maximum depth (of the top surface), etc. It returns depths in kilometers.
    The x, y coordinates are assumed to be in degrees. If using different model scaling, \
    this could be changed using kms_to_model
    
    This function will create layers that extend from the trench in the "y" direction. 
    It is assumed that the reference frame has been shifted so that relative plate motion is North-South, 
    i.e. the Y axis in python/numpy 
    
    Depths are considered positive
    """
    from scipy import interpolate
    #Create an interpolant for the location of the trench
    tpts = list(trench_shp.coords)
    px = [i[0] for i in tpts]
    py = [i[1] for i in tpts]
    f = interpolate.interp1d(px, py, kind='linear', axis=-1) #Linear interp just goes through the points exactly
    """
    Here we want a function that returns a function 
    The function returned is a function of x,y, f(x, y)
    """
    [i/kms_to_model for i in layer_depths]
    def create_surface(ldepth):
        rl = Rc - ldepth  #adjusted radius of curvature. 
        ang_rad = angle*(math.pi/180)
        slope = math.tan(ang_rad)  
        z_ec = math.sqrt((rl**2/ ((slope**2)+1))) #km
        y_ec = math.sqrt(rl**2 - z_ec**2)            #km
        Z_ec = abs(z_ec - Rc)  / 111.                 #model units, correct for the fact that z at theta = 0 is max, here we use Rc, not rl
        Y_ec = y_ec / 111.                            #model units          
        """
        Here we define the depth function f(x,y)
        """
        def surface(x,y):
            tpts = list(trench_shp.coords)
            px = [i[0] for i in tpts]
            py = [i[1] for i in tpts]
            f = interpolate.interp1d(px, py, kind='linear', axis=-1)
            px = x
            YaC = abs(y - float(f(px)))
            #print "YaC is %s." % YaC
            if YaC < Y_ec:              
                Ykm = YaC*111.
                #print "Rl is %s." % rl
                #print "Ykm is %s." % Ykm
                fd1 = math.sqrt(rl**2 - Ykm**2)  
                fd2 = abs(fd1 - rl)                  #correct for the fact that z at theta = 0 is max, here we use rl
                #print "fd2 is %s." % fd2
                fd3 = fd2 + ldepth                 #add the layer depth to get the total depth
                fdepth = fd3/111.
            if YaC > Y_ec:
                YaS = abs(y - float(f(px))) - Y_ec    #this guy is messing stuff up
                #YaS = y - float(f(px)) - Y_ec    #this guy is messing stuff up
                #print "Slope is %s." % slope
                fdepth = (YaS * slope) + Z_ec     #not sure about ldepth
                
            return fdepth*-111.
        return surface
        
    """
    Make the layers:
    """
    f_basket = []
    for dep in layer_depths:
        f = create_surface(dep)
        f_basket.append(f)
    return f_basket
    
    
    
    
def create_slab_region(trench_shp,direction = 1, Rc = 180, angle = 70, Zmax = 250, kms_to_model = 111.):
    """
    This function creates creates a polygon that covers the region
    where the slab will be interpolated, i.e. from the trench to the point
    defined by the parameters: 
    Radius of curvature, 
    Angle after which the slab is flat
    Zmax: the deepest part of the (top) slab
    """

    ang_rad = angle*(math.pi/180)
    slope = math.tan(ang_rad)  
    z_ec = math.sqrt((Rc**2/ ((slope**2)+1))) #km
    y_ec = math.sqrt(Rc**2 - z_ec**2)            #km
    Z_ec = abs(z_ec - Rc)  / kms_to_model                 #model units, correct for the fact that z at theta = 0 is max, here we use Rc, not rl
    Y_ec = y_ec / kms_to_model
    
    #work out northward distance covering the straight part of the slab
    tdep = Zmax - y_ec
    y2 = abs(tdep / slope)
    my2 = y2 / 111.
    
    tdis = my2 + Y_ec                      #model units
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
    
    
    
    
def map_flat_layers(position, dim, lon = [], lat = [], shapes=[], layers = [], upper_mantle_index=4, lower_mantle_index=7, lower_mantle_depth = 4, mantle_buffer = [10,-10,10,-10]):
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
                mat_counter = 0
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
    
    
    
def map_slab_layers(ii, mat_indx_var, position, dim, cutoff  = 250, lon = [], lat = [], shapes=[], layers = [], mantle_buffer = [10,-10,10,-10]):
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
    #pyGplates to be used in near-future
    else: 
        pt = Point(lon, lat)
        found_shape = 0
        mat_index = 0
        mat_counter = 0
        for i in range(0, len(shapes)):
            if i == 0:
                mat_counter = 0
            else: 
                mat_counter = mat_counter + ((len(layers[i-1])-1))
            mat_index = mat_counter
            found_layer = 0       
            if found_shape ==1:
                break
            elif shapes[i].contains(pt):
                found_shape == 1
                found_layer = 0 #reset after layer loop
                for j in range(0, len(layers[i][:-1])):   #This -2 thing would change depending on how many layers 
                    k = j+1
                    if found_layer == 1:
                        break
                    elif layers[i][j](lon,lat) >= depth > layers[i][k](lon,lat):
                        #print(layers[i][j](lon,lat))
                        mat_index = mat_index
                        found_layer = 1
                        break
                    elif k == len(layers[i]) -1:
                        mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0])
                        found_layer = 1
                        break
                    else:
                        mat_index = mat_index +1   
                break 
            else:
                mat_index = int(uw.swarms.tools.SwarmVariable_GetValueAt(mat_indx_var, ii)[0])
            
    return mat_index  