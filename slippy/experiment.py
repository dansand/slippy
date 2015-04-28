#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A very slapdash example of how you might build classes for parameter testing.
 - basic idea is you can create a base model, that is a, set of python instructions that
   when given appropriate paramaters, results in a valid (runnable) model
 - Experiments consist of series of models with varying parameters 
 - For each model in an experiment and unique list of input parameters needs to be stored,
    and a mode name that references the paramter(s) being modified in the experiment. 

:copyright:
    Dan Sandiford (sonderfjord@gmail.com), 2015
    Louis Moresi, 2015
:license:
    MIT
"""

###################
#Model setup
###################

#class Models(object):
#    def __init__(self):


#for i in ('apple', 'banana', 'carrot'):
#    fruitdict[i] = locals()[i]