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
       yield_prop=[50., 0.066667]
       
    .. rubric:: Example
    >>> test = materials.Materials(age = 120, crustThickness = 2)
    """
    def __init__(self, age = 80, crustThickness = 6, renorm= 2700, cooling_model="HSCM",
                 layers = [25, 50,75 ,100],
                 scales = [110000,1e20,10],
                 rock_prop=[2900, 3400,3500],
                 thermal_prop=[3.0e-5, 1e-6, 1300.0],
                 visc_prop=[240000, 5.0e-6, 5, 0.01],
                 yield_prop=[48, 0.66667, 5.0e-6]):
        
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
        ####Add some deeper values to layers:
        #This is a legacy of how the initial plotting stuff as done, needs to be fixed
        self.layers = layers
        self.layers.append(200)
        self.layers.append(100000)
        
        
    def write_properties_dict(self, out_name):
        import pickle
        dic = self.properties()
        n = out_name
        pickle.dump(dic, open( n, "wb" ))
                
    def properties(self, plot=False):
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
        
        
        ##################
        # Calculate values
        ##################
    
        for d in depth:
            if self.cooling_model == "HSCM":
                temp.append(utils.hsct(d*self.Kms,self.time,  self.Tint, self.Tsurf, self.kappa))
            elif self.cooling_model == "linear":
                temp.append(utils.linear_geotherm(d*self.Kms,self.layers[-3]*self.Kms,self.Tint, self.Tsurf))
                temp[index] = min(temp[index], self.Tint)
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
        #print LayerAverageDensity
        j = len(self.layers) - 2
        if len(self.scales) == 3:
            densityScale =  np.mean(LayerAverageDensity[0:j]) - self.hotMantleDensity
        else:
            densityScale = self.scales[3]
        print densityScale
        lds = np.mean(LayerAverageDensity[0:j]) #for the lithostatic scaling
        print lds
        #make sure no layer is less dense than the renorm value (assumes no thermal expansion in  top/crustal layer)
        #LayerAverageDensity = [max(i, self.renorm) for i in  LayerAverageDensity]
        #LayerAverageDensity = [((i - self.renorm) / densityScale) for i in LayerAverageDensity]
        LayerAverageDensity = [((i - self.renorm) / densityScale) for i in LayerAverageDensity]
        scaledMantleDensity = (self.hotMantleDensity - self.renorm) / densityScale
        #viscous time scale
        timeScale = self.viscosityScale/(densityScale* self.lengthScale * self.gravityScale)
        stokesstressScale = densityScale* self.gravityScale * self.lengthScale
        lithstressScale = lds* self.gravityScale * self.lengthScale/np.mean(LayerAverageDensity[0:j])
        #scales =[lengthScale,viscosityScale,gravityScale]
        self.scales.append(densityScale)
        self.scales.append(timeScale)
        self.scales.append(stokesstressScale)
        self.scales.append(lithstressScale)
        
        ##############
        # Make materials dict
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
        litho_dict["otherParams"] = [self.yield_prop[0], self.yield_prop[1],scaledMantleDensity]
        layer = self.layers
        litho_dict["layers"] = LayerTop[0:-1]
        scale = self.scales
        dodge = min(6, len(scale))
        litho_dict["scalingFactors"] = scale[0:dodge+1]
        


        
        ##############
        # Make plotting dict
        ##############  
        
        if plot == True:
            plot_dict = dict()
            plot_dict["LayerTop"] = LayerTop
            plot_dict["depth"] = depth
            plot_dict["temp"] = temp
            plot_dict["logViscosity"] = logViscosity
            plot_dict["logViscosity_d"] = logViscosity_d
            plot_dict["logViscosity_tr"] = logViscosity_tr
            plot_dict["density"] = density
            plot_dict["lithostaticPressure"] = lithostaticPressure
            plot_dict["yieldStrength"] = yieldStrength
            
            
            plot_dict["LayerLogAverageVisc"] = LayerLogAverageVisc
            plot_dict["LayerAverageDensity"] = LayerAverageDensity
            plot_dict["LayerAverageStrength"] = LayerAverageStrength
            

        # Return desired dictionary            
            return plot_dict
        else:
            return litho_dict
            
            
            
            
    def plot(self):
        """
        calls "properties" method for current instance, with plot flag set to True
        At this stage returns a default plot, adapted from Louis Moresi 
        """
        try:
            import matplotlib.pyplot as pypl
            from matplotlib.backends.backend_pdf import PdfPages
        except ImportError:
            print "The matplotlib tools are required to make plots with this script !"
            quit()
        
        #Gather any necessary class variables/attributes
        depth = self.layers
        LayerBase = self.layers
        barnum = len(depth) - 2
        
        
        
        #Get the necessary objects to plot, by calling self.properties
        pdict = self.properties(plot=True)
        
        #get the objects to plot 
        
        LayerTop = pdict["LayerTop"]
        depth = pdict["depth"]
        temp = pdict["temp"]
        logViscosity = pdict["logViscosity"]
        logViscosity_d = pdict["logViscosity_d"]
        logViscosity_tr = pdict["logViscosity_tr"]
        density = pdict["density"]
        lithostaticPressure = pdict["lithostaticPressure"]
        yieldStrength = pdict["yieldStrength"]
        
        LayerLogAverageVisc = pdict["LayerLogAverageVisc"]
        LayerAverageDensity = pdict["LayerAverageDensity"]
        LayerAverageStrength = pdict["LayerAverageStrength"]
        
        boxheight = [i - j for i,j in zip(LayerBase, LayerTop)]
        
        

        
        #Do plot
        figure1 = pypl.figure(figsize=(10,6))


        tempPlot = figure1.add_subplot(143)
        tempPlot.plot(temp,depth,  color='green', linestyle='solid', linewidth=2, label='')
        tempPlot.set_xlabel('Temperature ($^\circ$C)')
        tempPlot.set_xlim(0,1600)
        tempPlot.set_xticks([0,500,1000,1500])
        tempPlot.set_yticks([0,50,100,150, 200,250])
        tempPlot.grid(axis='y')
        tempPlot.invert_yaxis()
        pypl.setp( tempPlot.get_yticklabels(), visible=False)
               
        # now the viscosity plot

        viscPlot = tempPlot.twiny()
        viscPlot.plot(logViscosity_tr,depth, color='blue', linestyle='solid', linewidth=2, label='')
        viscPlot.plot(logViscosity,depth, color='blue', linestyle='dashed', linewidth=2, label='')
        viscPlot.plot(logViscosity_d,depth, color='#222222', linestyle='dashed', linewidth=2, label='')
        viscPlot.set_xlabel('log$_{10}$ relative viscosity')
        viscPlot.set_xlim(-2, 7)
        
        # now add viscosity-average bars
        viscPlot.barh(LayerTop[0:barnum],LayerLogAverageVisc[0:barnum], boxheight[0:barnum], left=0.0, alpha=0.5)  #Dont think this will cope with different width layers
        viscPlot.barh(LayerTop[barnum], LayerLogAverageVisc[barnum], boxheight[barnum], left=0.0, alpha=0.5, linestyle='solid', hatch='//')



        ## Add density plot

        densityPlot = figure1.add_subplot(142)
        densityPlot.set_xlabel('Density (kg/m$^3$)')
        
        densityPlot.set_xlim(2500,3500)
        densityPlot.set_xticks([2800,3000,3200,3400])
        densityPlot.set_yticks([0,50,100,150, 200,250])
        pypl.setp( densityPlot.get_yticklabels(), visible=False)

        densityPlot.grid(axis='y')
        densityPlot.plot(density,depth, color='red', linestyle='solid', linewidth=2, label='density (kg/m$^3$)')
        densityPlot.invert_yaxis()

        relDensityPlot = densityPlot.twiny()
        relDensityPlot.set_xlabel('Relative density')
        lhs = min(LayerAverageDensity[0:barnum])
        relDensityPlot.set_xlim(lhs-2, 20)    #make this a more resiliant
        relDensityPlot.axvline(x=10, ymin=0, ymax=1, linestyle='dotted', color='black')
        
        
        # now add relative density bars
        relDensityPlot.barh(LayerTop[0:barnum],LayerAverageDensity[0:barnum], boxheight[0:barnum], left=0.0, alpha=0.5, color='#FF2222')
        relDensityPlot.barh(LayerTop[barnum], LayerAverageDensity[barnum], boxheight[barnum], left=0.0, alpha=0.5, fill='false', color='#FF2222', hatch='//')


        ## Add strength/pressure plot
        pressurePlot = figure1.add_subplot(141)
        pressurePlot.set_xlabel('Pressure (Mpa)')
        # pressurePlot.set_xlim(2800,3300)
        pressurePlot.set_yticks([0,50,100,150, 200,250])
        pressurePlot.plot(lithostaticPressure, depth, color='DarkGreen', linestyle='dashed', linewidth=2, label='pressure (MPa)')
        pressurePlot.set_xticks([0,2500,5000,7500,10000])
        pressurePlot.grid(axis='y')
        pressurePlot.invert_yaxis()
        pressurePlot.set_ylabel('Depth (km)')

        pypl.setp( pressurePlot.get_yticklabels(), visible=True)
        
        
        ## Add the strength part 
        strengthPlot = pressurePlot.twiny()
        strengthPlot.set_xlabel('Yield strength (Mpa)')
        strengthPlot.plot(yieldStrength, depth, color='green', linestyle='solid', linewidth=2, label='Yield strength (MPa)')

        strengthPlot.barh(LayerTop[0:barnum],LayerAverageStrength[0:barnum], boxheight[0:barnum], left=0.0, alpha=0.5, color='green')
        # strengthPlot.barh(100.0, LayerAverageStrength[4], height=100.0, left=0.0, alpha=0.5, color='#EEFFEE', linestyle='solid', hatch='//')
        # Save this to 
        
        
        ##Add the combined strength plot
        #get the viscous part of the stress
        vel = 0.06 #velocity in m/y
        vels = vel / (365*24*60*60) #velocity in m per second
        visct = [(10**i)*self.viscosityScale for i in logViscosity_tr]
        visc = [(10**i)*self.viscosityScale for i in logViscosity]
        roc = 2.5e-11
        y = 50 * self.Kms
        exx = y*vels*roc
        visctStress =  [i*exx*1.0e-6 for i in visct]#viscous stress in MPa
        viscStress =  [i*exx*1.0e-6 for i in visc]#viscous stress in MPa
        
        
        cStrength = [min(i,j) for i,j in zip(viscStress, yieldStrength)]
        cStrengtht = [min(i,j) for i,j in zip(visctStress, yieldStrength)]
        
        
        #Do plot
        combstrengthPlot = figure1.add_subplot(144)
        combstrengthPlot.plot(cStrength,depth,  color='black', linestyle='solid', linewidth=2, label='')
        combstrengthPlot.plot(cStrengtht,depth, color='black', linestyle='dashed', linewidth=2, label='')
        combstrengthPlot.set_xticks([0,500, 1000,1500])
        combstrengthPlot.invert_yaxis()
        combstrengthPlot.set_ylim(150,0)
        
        

        pdf = PdfPages('viscosityDensityPlot-{:.0f}.pdf'.format(self.crustThickness))
           
        pypl.savefig(pdf, format='pdf') 
        pdf.close()
        
        out = viscStress[2425:2450]
        
        return LayerBase, LayerTop, out


        
             
        
            
            
        

        
            
               