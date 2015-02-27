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


def ocean_lith(age = 80, crustThickness = 10, cooling_model="HSCM",
    layers = [25,50,75 ,100,200, 1000000],
    scales = [111000,133,1e20,10,1.47e08, 6.79e11],
    rock_prop=[2900, 3400,3500],
    thermal_prop=[3.0e-5, 1e-6, 1300.0],
    visc_prop=[240000, 5.0e-6, 5, 0.01],
    yield_prop=[12.5, 0.066667, 5.0e-6]):
    '''
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
        scales = [lengthScale,densityScale,viscosityScale,stressScale] 
    -rock properties: densities (cold) for the rocks
        rock_prop=[basaltDensity (Oceanic crust - cold basalt), mantleDensity (cold), eclogiteDensity (cold)]
    -thermal properties: alpha - Thermal expansivity in /K, kappa - Thermal diffusivity m^2/s, upper mantle temp - degrees C]
        thermal_prop=[3.0e-5, 1e-6, 1300.0],
    -viscous properties: parameters for and Arrhenius-type relationship, Activation energy (J/mol),Activation volume, viscosityTruncation, viscosityMinimum
        visc_prop=[240000, 5.0e-6, 5, 0.01],
    -yield properties: parameters for the yielding rheology, cohesion, frictionCoeff
       yield_props=[12.5, 0.066667]    
    '''
    #Initial setup and constants
    TC2K = 273.0
    Tsurf = 0.0
    Tint = thermal_prop[2]
    kappa = thermal_prop[1]
    alpha = thermal_prop[0]
    Myrs = 3600*24*365.25*1.0e6     # Convert seconds to Myr
    Kms = 1000.0                    # Convert m to km
    R0 = 8.314                      # Gas constant
    gravAcc = 9.8
    E0 = visc_prop[0]
    crustDensity = rock_prop[0]
    mantleDensity = rock_prop[1]
    referenceDensityDifference = mantleDensity * (Tint-Tsurf) * alpha
    hotMantleDensity = mantleDensity - referenceDensityDifference 
    V0 =  mantleDensity * gravAcc * visc_prop[1] * 1.0 
    gravAcc = 9.8
    ap = mantleDensity*gravAcc
    cohesion = yield_prop[0]
    frictionCoeff = yield_prop[1]
    
    # rheology & temperature profile functions
    def viscosity(T,D):
        return math.exp((E0 + V0*D)/(R0*T))
    
    if cooling_model=="HSCM":
        "time is in years"
        def hsct(depth,time): 
            T = Tsurf +(Tint-Tsurf) * math.erf(depth/(2*math.sqrt(time*kappa)))
            return T
    
    #Other parameters
    time = age*Myrs
    viscosityTruncation = visc_prop[2]
    viscosityMinimum = visc_prop[3]
    stressScale = scales[0]* scales[1] * scales[3]
    referenceViscosity = viscosity(Tint+TC2K,100*Kms);
    #Create layers to iterate over
    # Layer definitions (Include a final layer to accumulate any out-of-range values used for plotting )
    depth = np.linspace(0.0, 250, num=5000) 
    temp = []
    density = []
    lithostaticPressure = []
    yieldStrength = []
    logViscosity = []
    logViscosity_tr = []
    logViscosity_d = []


    LayerBase  = layers 
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
    
    ##############
    # Calculate values
    ##############
    
    for d in depth:
        temp.append(hsct(d*Kms,time))
        density.append(DensityLayer[densitylayer] * (1 - (temp[index] - Tsurf) * alpha ))
        logViscosity.append(math.log(max(viscosity(temp[index]+TC2K,d*Kms),viscosityMinimum)/referenceViscosity,10.0))
        logViscosity_d.append(math.log(viscosity(Tint+TC2K,d*Kms)/referenceViscosity,10.0))
        logViscosity_tr.append(min(logViscosity[index],viscosityTruncation))
        
        
        
        if (index > 0):
        # lithospheric pressure in MPa
            lithostaticPressure.append(1.0e-6 * gravAcc  * (depth[index] - depth[index-1]) * Kms * density[index] + lithostaticPressure[index-1])
            yieldStrength.append(cohesion + frictionCoeff * lithostaticPressure[index])
            
        LayerAverageTemp[layer] = LayerAverageTemp[layer] + temp[index]
        LayerAverageDensity[layer] = LayerAverageDensity[layer] + density[index]
        LayerAverageVisc[layer] = LayerAverageVisc[layer] + viscosity(temp[index]+TC2K,d*Kms)/referenceViscosity
        LayerAverageStrength[layer] = LayerAverageStrength[layer] + lithostaticPressure[index] * frictionCoeff

        samples = samples + 1
        
        if d > LayerBase[layer]:
        #print '{} points in layer {} to depth {}'.format(samples,layer,d )

            LayerAverageTemp[layer] = LayerAverageTemp[layer] / samples
            LayerAverageDensity[layer] = ((LayerAverageDensity[layer] / samples) - hotMantleDensity)/ referenceDensityDifference
            LayerAverageVisc[layer] = LayerAverageVisc[layer] / samples
            LayerLogAverageVisc[layer] = min(viscosityTruncation,math.log(LayerAverageVisc[layer],10.0))
            LayerAverageStrength[layer] = cohesion + LayerAverageStrength[layer] / samples  # in MPa

            layer = layer + 1
            samples = 0	

        if d > DensityLayerBase[densitylayer]:
            densitylayer = densitylayer + 1	       

        
        index = index + 1
    
    ##############
    # Write layer averages to dictionary
    ##############
    
    litho_dict = dict()
    for i in range(0, len(layers)-2):
        layer1strength = cohesion * 1.0e6 + frictionCoeff * mantleDensity * gravAcc * (LayerTop[0] + LayerBase[0]) * 0.5 * Kms
        name = "layer" + str(i)
        litho_dict[name] = dict()      
        litho_dict[name]["AverageTemp"] = LayerAverageTemp[i]
        litho_dict[name]["AverageVisc"] = min(LayerAverageVisc[i], 10**(viscosityTruncation))
        litho_dict[name]["AverageDensity"] = LayerAverageDensity[i]
        litho_dict[name]["AverageStrength"] = layer1strength/stressScale
    litho_dict["scalingFactors"] = scales
    litho_dict["otherParams"] = [yield_prop[0], yield_prop[1]]
    return litho_dict


