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
    trench = uplate.intersection(lplate)
    #trench1 = translate(trench, xoff=0.0, yoff=my, zoff=0.0)
    trench2 = translate(trench, xoff=0.0, yoff=tdis, zoff=0.0)
    tpts = list(trench.coords)
    #Create shapefiles for the curved and straight slab part
    tpts2 = list(trench2.coords)
    tpts2.reverse()
    temp_crds = tpts + tpts2
    slab_region = Polygon(temp_crds)
    return(slab_region)