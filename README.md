# Scintillation Tutorial

The idea of this tutorial is to create a jupyter notebook to explore and visualize the basic mechanisms and observables of the scintillation of a compact radio source.

## Setup

Python3 with these libraries installed is required:
-  notebook
-  numpy
-  matplotlib

First, download *scintillation_tutorial.py* from this repository.
Then start the jupyter notebook in the directory that contains this file:
~~~
cd PATH/TO/DIRECTORY
jupyter notebook
~~~

Within the notebook, we start by importing the libraries used:
~~~
import numpy as np
import matplotlib.pyplot as plt
from scintillation_tutorial import ScatteredSignal
~~~

Next, we define the numerical values of some units conversions:
~~~
#constants
au = 149597870700. #m
pc = 648000./np.pi*au #m
day = 24*3600 #s
year = 365.2425*day #s
degrees = np.pi/180. #rad
mas = degrees/1000./3600. #rad
kHz = 1.0e+3 #Hz
MHz = 1.0e+6 #Hz
mus = 1.0-6 #s
minute = 60. #s
kms = 1000. #m/s
~~~
Throughout this tutorial, the values of all variables with physical dimensions will be set in units of meter, second, radian, and hertz to avoid errors.

## Single scattered rays

We start by creating a very simple scattering geometry. First, we need to define the distance of the source. To get an intuition, use the distance of a pulsar of FRB that you know. The example shown here follows numbers that are reasonable for PSR B1508+55.
The following lines create an instance of the *ScatteredSignal* class with your chosen distance:
~~~
Sc = ScatteredSignal(D = 2100.*pc)
~~~
Now we define a scattered path by introducing a point in space that is passed by an additional ray from the source to the observer. Again, use numbers of your choice:
~~~
Sc.addPoint(D=125.*pc,x=0.1*au)
~~~
This class method adds a new scattering point. If there are already points, they will not be deleted.
In this example, the pulsar is 2100 parsecs away. At a distance of 125 pc from Earth, an unspecified object that is offset from the direct line of sight by 0.1 astronomical units scatters the radio waves.
