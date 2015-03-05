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