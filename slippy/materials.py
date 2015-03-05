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



class Materials(object):
    """
    This function creates a dictionary with model inputs for undworld,
    Including, density, viscosity, temperature, yeilding.
    At the moment, the first object of this Class created needs to be the slab-
    This sets the scaling for the models, an other components follow. 
    ...
    -age: age of lithosphere (Myr)
    -crustThickness: km
    -renorm: this should be the density of the least dense part of your model
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
       
    .. rubric:: Example
    >>> test = materials.Materials(age = 120, crustThickness = 2)
    """
    def __init__(self, age = 80, crustThickness = 6, renorm= 2700, cooling_model="HSCM",
                 layers = [25,50,75 ,100,200, 1000000],
                 scales = [110000,1e20,10],
                 rock_prop=[2900, 3400,3500],
                 thermal_prop=[3.0e-5, 1e-6, 1300.0],
                 visc_prop=[240000, 5.0e-6, 5, 0.01],
                 yield_prop=[48, 0.066667, 5.0e-6]):
        
        if cooling_model == "HSCM":
            self.cooling_model = "HSCM"
        elif cooling_model =="linear":
            self.cooling_model = "linear"
        else:
            raise ValueError("Unknown cooling model '%s'." %
                             cooling_model)                           
        self.age = age
        self.crustThickness = crustThickness
        self.renorm = renorm
        self.cooling_model = cooling_model
        self.layers = layers
        self.scales = scales
        self.thermal_prop = thermal_prop
        self.visc_prop = visc_prop
        self.yield_prop = yield_prop
        self.rock_prop = rock_prop
        #Model scales
        self.lengthScale = scales[0]
        self.viscosityScale = scales[1]
        self.gravityScale = scales[2]
        #Physical Constants
        self.TC2K = 273.0
        self.Kms = 1000.0
        self.gravAcc = 9.8
        self.R0 = 8.314
        #Other vars
        self.Tsurf = 0.0
        self.Tint = thermal_prop[2]
        self.kappa = thermal_prop[1]
        self.alpha = thermal_prop[0]
        self.cohesion = yield_prop[0]
        self.frictionCoeff = yield_prop[1]
        self.E0 = visc_prop[0]
        self.crustDensity = rock_prop[0]
        self.mantleDensity = rock_prop[1]
        self.V0 =  self.mantleDensity * self.gravAcc * self.visc_prop[1] * 1.0 
        self.Tsurf = 0.0
        self.Tint = thermal_prop[2]
        self.referenceDensityDifference = self.mantleDensity * (self.Tint-self.Tsurf) * self.alpha
        self.hotMantleDensity = self.mantleDensity - self.referenceDensityDifference
        self.Myrs = 3600*24*365.25*1.0e6
        self.time = self.age*self.Myrs
        self.viscosityTruncation = visc_prop[2]
        self.viscosityMinimum = visc_prop[3]
        self.referenceViscosity = utils.viscosity(self.Tint+self.TC2K,100*self.Kms, self.E0, self.V0, self.R0);
        
    def dan_func(self, param):
     
        self.viscosityTruncation = self.visc_prop[2]
        crazy_val2 = 6.
        some_value = self.Tint
        out = utils.viscosity(self.Tint+self.TC2K,100*self.Kms, self.E0, self.V0, self.R0)
        out = out + self.viscosityTruncation + crazy_val2
        return out
        
    def write_properties_dict(self, out_name):
        import pickle
        dic = self.properties()
        n = out_name
        pickle.dump(dic, open( n, "wb" ) )
        
        
        
    def properties(self):
        # Layer definitions (Include a final layer to accumulate any out-of-range values used for plotting )
        depth = np.linspace(0.0, 250, num=5000) 
        temp = []
        density = []
        lithostaticPressure = []
        yieldStrength = []
        logViscosity = []
        logViscosity_tr = []
        logViscosity_d = []
        LayerBase  = self.layers 
        LayerTop = [0] + LayerBase[0:len(LayerBase)-1]
        LayerAverageTemp = [0,0,0,0,0,0]
        LayerAverageDensity = [0,0,0,0,0,0]
        LayerAverageStrength = [0,0,0,0,0,0]
        LayerAverageVisc = [0,0,0,0,0,0]
        LayerLogAverageVisc = [0,0,0,0,0,0]
        DensityLayerBase = [self.crustThickness,1000000]
        DensityLayer = [self.crustDensity,self.mantleDensity]
        index = 0
        layer = 0
        samples = 0
        densitylayer = 0
        lithostaticPressure.append(0.0)
        yieldStrength.append(self.cohesion)
        ##############
        # Calculate values
        ##############
    
        for d in depth:
            if self.cooling_model == "HSCM":
                temp.append(utils.hsct(d*self.Kms,self.time,  self.Tint, self.Tsurf, self.kappa))
            elif self.cooling_model == "linear":
                temp.append(utils.linear_geotherm(d*self.Kms,self.layers[-3]*self.Kms,self.Tint, self.Tsurf))
            density.append(DensityLayer[densitylayer] * (1 - (temp[index] - self.Tsurf) * self.alpha ))
            logViscosity.append(math.log(max(utils.viscosity(temp[index]+self.TC2K,d*self.Kms, self.E0, self.V0, self.R0),self.viscosityMinimum)/self.referenceViscosity,10.0))
            logViscosity_d.append(math.log(utils.viscosity(self.Tint+self.TC2K,d*self.Kms, self.E0, self.V0, self.R0)/self.referenceViscosity,10.0))
            logViscosity_tr.append(min(logViscosity[index],self.viscosityTruncation))
            
            if (index > 0):
            # lithospheric pressure in MPa
                lithostaticPressure.append(1.0e-6 * self.gravAcc  * (depth[index] - depth[index-1]) * self.Kms * density[index] + lithostaticPressure[index-1])
                yieldStrength.append(self.cohesion + self.frictionCoeff * lithostaticPressure[index])
            
            LayerAverageTemp[layer] = LayerAverageTemp[layer] + temp[index]
            LayerAverageDensity[layer] = LayerAverageDensity[layer] + density[index]
            LayerAverageVisc[layer] = LayerAverageVisc[layer] + utils.viscosity(temp[index]+self.TC2K,d*self.Kms, self.E0, self.V0, self.R0)/self.referenceViscosity
            LayerAverageStrength[layer] = LayerAverageStrength[layer] + lithostaticPressure[index] * self.frictionCoeff

            samples = samples + 1
        
            if d > LayerBase[layer]:
            #print '{} points in layer {} to depth {}'.format(samples,layer,d )

                LayerAverageTemp[layer] = LayerAverageTemp[layer] / samples
                LayerAverageDensity[layer] = (LayerAverageDensity[layer] / samples)
                LayerAverageVisc[layer] = LayerAverageVisc[layer] / samples
                LayerLogAverageVisc[layer] = min(self.viscosityTruncation,math.log(LayerAverageVisc[layer],10.0))
                LayerAverageStrength[layer] = self.cohesion + LayerAverageStrength[layer] / samples  # in MPa

                layer = layer + 1
                samples = 0	

            if d > DensityLayerBase[densitylayer]:
                densitylayer = densitylayer + 1	
        
            index = index + 1
            
            
        ##############
        # Rescaling
        ##############
        j = len(self.layers) - 2
        if len(self.scales) == 3:
            densityScale =  np.mean(LayerAverageDensity[0:j]) - self.hotMantleDensity
        else:
            densityScale = self.scales[3]
        LayerAverageDensity = [(i - self.renorm) / densityScale for i in LayerAverageDensity]
        scaledMantleDensity = (self.hotMantleDensity - self.renorm)/densityScale
        timeScale = self.viscosityScale/(densityScale* self.lengthScale * self.gravityScale)
        stokesstressScale = densityScale* self.gravityScale * self.lengthScale
        lithstressScale = densityScale* self.gravityScale * self.lengthScale * (np.mean(LayerAverageDensity[0:4]))
        #scales =[lengthScale,viscosityScale,gravityScale]
        self.scales.append(densityScale)
        self.scales.append(timeScale)
        self.scales.append(stokesstressScale)
        self.scales.append(lithstressScale)
        
        ##############
        # Rescaling
        ##############
        
        litho_dict = dict()
        for i in range(0, len(self.layers)-2):
            #layer1strength = cohesion * 1.0e6 + frictionCoeff * mantleDensity * gravAcc * (LayerTop[0] + LayerBase[0]) * 0.5 * Kms
            name = "layer" + str(i)
            litho_dict[name] = dict()      
            litho_dict[name]["AverageTemp"] = LayerAverageTemp[i]
            litho_dict[name]["AverageVisc"] = min(LayerAverageVisc[i], 10**(self.viscosityTruncation))
            litho_dict[name]["AverageDensity"] = LayerAverageDensity[i]
            litho_dict[name]["AverageStrength"] = (LayerAverageStrength[i] * 1.0e6)/lithstressScale
        litho_dict["scalingFactors"] = self.scales
        litho_dict["otherParams"] = [self.yield_prop[0], self.yield_prop[1],scaledMantleDensity]
        return litho_dict    

        
            
               