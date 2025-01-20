# ScintillationTutorial

The idea of this tutorial is to create a jupyter notebook to explore and visualize the basic mechanisms and observables of the scintillation of a compact radio source.

## Setup

Python3 with these libraries installed is required:
-  notebook
-  numpy
-  matplotlib

First, download *scintillation_tutorial.py* from this repository.
Then start the jupyter notebook in the directory that contains this file:
> cd PATH/TO/DIRECTORY
> jupyter notebook

Within the notebook, we start by importing the libraries used:
> import numpy as np
> import matplotlib.pyplot as plt
> from scintillation_tutorial import ScatteredSignal

Next, we define the numerical values of some constants and units:
> #constants
> au = 149597870700. #m
> pc = 648000./np.pi*au #m
> day = 24*3600 #s
> year = 365.2425*day #s
> degrees = np.pi/180. #rad
> mas = degrees/1000./3600. #rad
> kHz = 1.0e+3 #Hz
> MHz = 1.0e+6 #Hz
> mus = 1.0-6 #s
> minute = 60. #s
> kms = 1000. #km/s



