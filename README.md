# slippy

Slippy is a small python package to help run geodynamic models using the Underworld framework. Primarily, it will help integrate Pygplates (currently unreleased) and Underworld2 (Python version of Underworld). It is hoped that in the future the only dependencies will be pygplates and underworld (python package). For the time being, it will rely on a couple of other packages, particularly Shapely, to fill out some topology requirements.

The main branch of this repository relates to a deprecated version of Underworld2; the newInterface branch corresponds to a branch of the same name in the Underworld project where all development is currently focussed. This will become the main branch of Slippy.

to install slippy:


    git clone https://github.com/dansand/slippy

    cd slippy

    pip install -v -e .
