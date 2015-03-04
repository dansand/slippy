#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Materials class of slippy.

:copyright:
    Dan Sandiford (sonderfjord@gmail.com), 2015
    Louis Moresi, 2015
:license:
    MIT
"""


import math
import numpy as np

from slippy import utils



class CreateMaterials(object):
    """
    This function creates a dictionary with model inputs for undworld,
    Including, density, viscosity, temperature, yeilding
    ...
    -age: age of lithosphere (Myr)
    -crustThickness: km
    - cooling model: The cooling model is used to derive densities.
        cooling_model="HSCM", "PLATE"
    -layers: The user can specify any number of layers, these are the depths (km) of the bottom of each layer
        layers = [25,50,75 ,100,200, 1000000],
    -scales: scaling values for the system
        scales = [lengthScale,viscosityScale,gravityScale] 
    -rock properties: densities (cold) for the rocks
        rock_prop=[basaltDensity (Oceanic crust - cold basalt), mantleDensity (cold), eclogiteDensity (cold)]
    -thermal properties: alpha - Thermal expansivity in /K, kappa - Thermal diffusivity m^2/s, upper mantle temp - degrees C]
        thermal_prop=[3.0e-5, 1e-6, 1300.0],
    -viscous properties: parameters for and Arrhenius-type relationship, Activation energy (J/mol),Activation volume, viscosityTruncation, viscosityMinimum
        visc_prop=[240000, 5.0e-6, 5, 0.01],
    -yield properties: parameters for the yielding rheology, cohesion, frictionCoeff
       yield_props=[50., 0.066667]    
    """
    def __init__(self, age = 80, crustThickness = 6, cooling_model="HSCM",
                 layers = [25,50,75 ,100,200, 1000000],
                 scales = [lengthScale,viscosityScale,gravityScale],
                 rock_prop=[2900, 3400,3500],
                 thermal_prop=[3.0e-5, 1e-6, 1300.0],
                 visc_prop=[240000, 5.0e-6, 5, 0.01],
                 yield_prop=[48, 0.066667, 5.0e-6]):
        self.age = age
        self.crustThickness = crustThickness
        self.cooling_model = cooling_model
        self.layers = layers
        self.scales = scales
        self.thermal_prop = thermal_prop
        self.visc_prop = visc_prop
        self.yield_prop = yield_prop
        #Initial setup and constants
        self.TC2K = 273.0
        self.Tsurf = 0.0
        self.Tint = self.thermal_prop[2]
        self.kappa = self.thermal_prop[1]
        self.alpha = self.thermal_prop[0]
        self.Myrs = 3600*24*365.25*1.0e6     # Convert seconds to Myr
        self.Kms = 1000.0                    # Convert m to km
        self.gravAcc = 9.8
        self.E0 = visc_prop[0]
        self.crustDensity = rock_prop[0]
        self.mantleDensity = rock_prop[1]
        self.referenceDensityDifference = mantleDensity * (Tint-Tsurf) * alpha
        self.hotMantleDensity = mantleDensity - referenceDensityDifference 
        self.V0 =  mantleDensity * gravAcc * visc_prop[1] * 1.0 
        self.R0 = 8.314                      # Gas constant
        self.gravAcc = 9.8
        self.ap = mantleDensity*gravAcc
        self.cohesion = yield_prop[0]
        self.frictionCoeff = yield_prop[1]  
    #Other parameters
        self.time = self.age*self.Myrs
        self.viscosityTruncation = visc_prop[2]
        self.viscosityMinimum = visc_prop[3]
        self.referenceViscosity = viscosity(Tint+TC2K,100*Kms);
       #Create layers to iterate over
       # Layer definitions (Include a final layer to accumulate any out-of-range values used for plotting )
       
        def test_function(self):
            depth = np.linspace(0.0, 250, num=5000) 
            temp = []
            density = []
            lithostaticPressure = []
            yieldStrength = []
            logViscosity = []
            logViscosity_tr = []
            logViscosity_d = []
            #More listst
            LayerBase  = self.layers 
            LayerTop = [0] + LayerBase[0:len(LayerBase)-1]
            LayerAverageTemp = [0,0,0,0,0,0]
            LayerAverageDensity = [0,0,0,0,0,0]
            LayerAverageStrength = [0,0,0,0,0,0]
            LayerAverageVisc = [0,0,0,0,0,0]
            LayerLogAverageVisc = [0,0,0,0,0,0]
            DensityLayerBase = [crustThickness,1000000]
            DensityLayer = [crustDensity,mantleDensity]
    
            index = 0
            layer = 0
            samples = 0
            densitylayer = 0
            lithostaticPressure.append(0.0)
            yieldStrength.append(cohesion)
            
            return LayerBase
  