# slippy

Slippy is a small python package to help steer geodynamic models using the Underworld framework. At the moment it relies on the python package Shapely, but will likely switch to using pygplates for the necessary functions. 

At this stage, slippy allows you to create models with planform complexity: basically you set up a set of polygons, then for each polygon you set up a series of layers and material properties. 

The workflow currently only supports visco-plastic problems, i.e. the heat equation is not solved. Temperature dependent material properties are estimated using inferred temperature profiles (like halfspace cooling solutions).

The map_layers functions then loop through shapes, and layers and assign the appropriate material properties to the particle swarm. 

One thing to watch out for at the moment is that the oceanic lithosphere (i.e. the active part of the system) should be the first material in the shapes and layers list. This is because some system scalings are set using the properties of the slab, then impoosed on other parts of the model.

to install slippy:


    git clone https://github.com/dansand/slippy
    
    cd slippy
    
    pip install 
   

