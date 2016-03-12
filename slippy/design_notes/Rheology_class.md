## Rheology

Need a good convection for notation of dimensional and non-dimensional coordinates:

In the future, you will be able to import a rheology by using the Rheology class.

The Basic rheology class will consist of:

*  dimensional units / dimensionless units
* Scaling factors

Such that, taking either dimensionless or dimensional values from a paper, we can create the required inputs. 

`print Rheology.avalable_models()`

Individual Rheologies will inherit from the Rheology class, taking things like scaling factors (as these are fairly standard),

`model_rheology = Rheology.Tosi2015(temperatureField, DepthField, secInvStrainRate)`

A particular rheology will have a number of materials (hopefully one day just 1 for everything)

``` model_rheology.avalable_materials()
model_rheology.available_materials()
	$ uppermantle
```



``` model_rheology.available_models()	$ mantle_lithoshere
modelX_rheology.available_materials()
	$ lowermantle
    $ uppermantle
    $ crust
    $ air
```

We may have to have a `material class, so we can do:`

`lowermantle = model_rheology.lowermantle()`

`lowermantlefunction = lowermantle.createfuntion()`

If you actually want to grab the parameter values you could do:

`lowermantle.dimonsional.cohesion()`

`lowermantle.nondimensional.cohesion()`

Where these read from dictionaries, i.e, in the class definition

`lowermantle.nondimensional['cohesion']`

Or something like that...

The dimensional and non dimensional dictionaries may not map 1-1, for instance, the density and gravity... etc paramaters in the dimensional dictionare effectively collaps into the Rayleigh number in the nondimensiona dict. 









The rheology class will then produce uw.functions for different rheologies…These should be available as 

`model_rheology.uppermantle()`

`model_rheology.lithosphere()`

Other methods.

`model_rheology.print_reference()`

​	Bibtex of paper

`model_rheology.print_dimensional_parameters()`

`model_rheology.print_dimensionless_parameters()`

### Syntax

How should we denote dimensionless vs dimensional variables?

model_rheology.dimensional.eta = 

model_rheology.nondimensional.eta = 