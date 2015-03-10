#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Shapes class of slippy.

:copyright:
    Dan Sandiford (sonderfjord@gmail.com), 2015
    Louis Moresi, 2015
:license:
    MIT
"""
import math
import numpy as np
import underworld as uw



def passive_tracer_carpet(minX=minX, maxX=maxX, minY=minY, maxY=maxY,minZ=minZ,maxZ=maxZ, numx=10, numy=10, swidth=0.5):
    '''
    Create a passive tracer swarm in a cross-hatch shape
    Default args assume that domain boundaries are e.g. minX, maxX; 
    Replace those if necessary
    '''
    xshapes_names = []
    yshapes_names = []
    #Create the Y-parallel shapes
    if numx != 0:
        #create list of coordinate values (corners) for strips parallel to Y axis:
        llist = []
        lrlist = []
        dx = (abs(maxX) + abs(minX))/(numx + 1)
        #print dx
        lc = minX + (dx - 0.5*swidth)
        rc = minX + (dx + 0.5*swidth)
        #print lc, rc
        for i in range(1, numx+1):
            llist.append(lc)
            lrlist.append(rc)
            lc = lc + dx
            rc = rc + dx
            name = 'yboxShape' + str(i)
            yshapes_names.append(name)
            shape.setup.boxCreate(componentName=name, startList=[llist[i-1],minY,minZ], endList=[lrlist[i-1],maxY,maxZ])
        #Create the X-parallel shapes
        #print llist, lrlist
    if numy != 0:
        #create list of coordinate values (corners) for strips parallel to X axis:
        blist = []
        tlist = []
        dy = (abs(maxY) + abs(minY))/(numy + 1)
        bc = minY + (dy - 0.5*swidth)
        tc = minY + (dy + 0.5*swidth)
        for i in range(1, numy+1):
            blist.append(bc)
            tlist.append(tc)
            bc = bc + dy
            tc = tc + dy
            name = 'xboxShape' + str(i)
            xshapes_names.append(name)
            shape.setup.boxCreate(componentName=name, startList=[minX,blist[i-1],minZ], endList=[maxX,tlist[i-1],maxZ])
    out = dict([('Type', 'Union'), ('shapes', yshapes_names + xshapes_names)])
    return out  