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

We start by creating a very simple scattering geometry. First, we need to define the distance of the source. To get an intuition, use the distance of a pulsar of FRB that you know. In the case of FRBs, all distances are angular diameter distances. The example shown here follows numbers that are reasonable for PSR B1508+55.
The following lines create an instance of the *ScatteredSignal* class with your chosen distance:
~~~
Sc = ScatteredSignal(D = 2100.*pc)
~~~
Now we define a scattered path by introducing a point in space that is passed by an additional ray from the source to the observer. Again, use numbers of your choice:
~~~
Sc.addPoint(D=125.*pc,x=0.1*au)
~~~
This class method adds a new scattering point. If there are already points, they will not be deleted.\
In this example, the pulsar is 2100 parsecs away. At a distance of 125 pc from Earth, an unspecified object that is offset from the direct line of sight by 0.1 astronomical units scatters the radio waves.
You can plot this scattering geometry with
~~~
Sc.plot_rays()
~~~
This plot shows the paths the radio waves are taking from the source at the right to the observer on the left. Of course, these paths are not directly observable. For this reason, the main work within the field of scintillometry is to reconstruct them from the data we get.\
If we have instruments with sufficient angular resolution, we can observe that the scattered rays arrive from different angular positions on the sky. Such a plot of the sky can be created with
~~~
Sc.plot_images()
~~~
Note that this plot assumes that the source is a point. For an extended source each of the scatterd paths would look like the source and potentially overlap to the point where they cannot be separated. In scintillometry jargon, we call each scattered path an *image* for this reason since each path creates another image of the source on the sky. Much of the terminology is shared with strong gravitational lensing. But keep in mind that in the case of scintillation, each of these images is usually assumed to be just a point because the source is compact. If there is a large cloud of images, the scattered source may appear like a *scattering disk*.\
To enter an image that is offset in both dimensions on the sky, use
~~~
Sc.addPoint(D=425.*pc,x=-0.1*au,y=0.5*au)
~~~
