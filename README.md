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
mus = 1.0e-6 #s
minute = 60. #s
kms = 1000. #m/s
~~~
Throughout this tutorial, the values of all variables with physical dimensions will be set in units of meter, second, radian, and hertz to avoid errors.

## Single scattered rays

We start by creating a very simple scattering geometry. First, we need to define the distance of the source. To get an intuition, use the distance of a pulsar, FRB, or AGN that you know. In the case of extragalactic objects, all distances are angular diameter distances. The example shown here follows numbers that are reasonable for PSR B1508+55.
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
This plot shows the paths the radio waves are taking from the source at the right to the observer on the left. Of course, these paths are not directly observable. For this reason, the main work within the field of scintillometry is to reconstruct them from the data we get.

### Scattering disk

If we have instruments with sufficient angular resolution, we can observe that the scattered rays arrive from different angular positions on the sky. Such a plot of the sky can be created with
~~~
Sc.plot_images()
~~~
Note that this plot assumes that the source is a point. For an extended source each of the scatterd paths would look like the source and potentially overlap to the point where they cannot be separated. In scintillometry jargon, we call each scattered path an *image* for this reason since each path creates another image of the source on the sky. Much of the terminology is shared with strong gravitational lensing. But keep in mind that in the case of scintillation, each of these images is usually assumed to be just a point because the source is compact. If there is a large cloud of images, the scattered source may appear like a *scattering disk*.\
To enter an image that is offset in both dimensions on the sky, use
~~~
Sc.addPoint(D=425.*pc,x=-0.1*au,y=0.5*au)
~~~
Compare the size of your scattering disk to the $\lambda / D$ of a large dish or the baselines possible in ground-based VLBI to see if you could observe it.

### Scattering delay

Pulsed sources like pulsars allow for the measurement of the time of arrival of radiation. Hence, parts of the signal that are delayed can be separated. Signals are delayed because of the changed speed of propagation within a medium (dispersive delay) and because of different lengths of their paths between source and observer (geometric delay). The dispersion delay of the bulk of the system is what comprises the dispersion measure (DM). While each path could have a slightly different DM, this has, interestingly, not been observed to be significant. The geometric delay between the paths is much larger, so we will neglect the dispersive delay.\
The delay is obtained by dividing the difference between the length of the scattered path and the direct path by the speed of light. In the small-angle approximation the triangular shape leads to the delay being proportional to the square of the offset of the scatterer from the direct line of sight.\
You can look at the delays of your system with
~~~
Sc.plot_delays()
~~~
Again, the intrinsic signal is assumed to be point-like. In reality, a pulse or burst will have a nonzero length and structure that gets convolved with the delays, similar to the case of the images above.

### Scintillation

Scintillation is an interference effect. The waves traveling along different scattered paths arrive also with phase delays with respect to each other. Stable interference patterns are only possible for coherent radiation. The phase evolution of waves depends on the radio frequency such that a particular spectral pattern is obtained.\
To simulate scintillation, we first need to define a frequency band, for example:
~~~
nu = np.linspace(1200.*MHz,1600.*MHz,num=200)
~~~
Since the computation can be longer for large numbers of images and is not very optimized here, it makes sense to separate computation and plotting. Scintillation is computed with
~~~
I = Sc.compute_scintillation(nu)
~~~
and plotted with
~~~
Sc.plot_scintillation(nu,I)
~~~
Each pair of interfering paths introduces another sinusoidal pattern to the scintillation. As a result, a low number of scattered paths leads to a clearly periodic pattern while a patchy random pattern without clear periodicity is indicative of a high number of scattered paths.\

In the first case, scintillation is much easier to study in Fourier space, while the second case is intrinsically more difficult to study. In fact, the vast majority of modern pulsar scintillometry happens solely in Fourier space, which is why it is important to familiarize yourself with it.\
The power spectrum is a good representation of the strength of each conjugate frequency component. It is computed as the squared modulus of the Fourier transform:
~~~
power = np.abs(np.fft.fft(I))**2
~~~
The values of the conjugate frequencies are obtained by
~~~
N_nu = len(nu)
dnu = nu[1]-nu[0]
f_nu = np.fft.fftfreq(N_nu,dnu)
~~~
These coordinates are not sorted from smallest to highest value by definition. This needs to be changed:
~~~
power = np.fft.fftshift(power)
f_nu = np.fft.fftshift(f_nu)
~~~
Now, we can plot the power spectrum:
~~~
figure = plt.figure(figsize=(16,9))
ax = figure.add_subplot(1,1,1)
ax.plot(f_nu*MHz,power/np.max(power),color="black",marker="",linestyle='-')
ax.set_xlabel(r"Delay $f_{\nu}=\tau$ [1/MHz=µs]")
ax.set_ylabel(r"Power $\vert {\rm FT}(I)\vert^2$ (max=1)")
~~~
The range of values and the resolution within a power spectrum are a direct consequence of the value range and channel width of the observed frequencies. They are inverse to each other:
~~~
BW = nu[-1]-nu[0]
print("Number of channels: {0}\nChannel width: {1:.2e} MHz\tinverted: {2:.2e} µs\nBandwidth: {3:.2e} MHz\t\tinverted: {4:.2e} µs".format(N_nu,dnu/MHz,1/(dnu/MHz),BW/MHz,1/(BW/MHz)))
~~~
The power spectrum is affected by many Fourier artifacts for which different mitigation techniques have been developed. A simple one is to subtract the mean before performing the FFT, which removes the central bright peak:
~~~
power = np.abs(np.fft.fft(I-np.mean(I)))**2
~~~
The coordinates of spikes of power can be predicted from the image positions defined earlier and represent a differential version of their delays. That's why the Fourier conjugate of frequency is also called delay $\tau$. We can compare the theoretical positions of power:
~~~
Sc.plot_differential_delays(y_max=0.1)
~~~
The naming can be confusing because this is the delay difference between all pairs of images instead of between one image and the direct line of sight. Also, we look at the phase delay here, while before we were looking at the group delay. Both happen to be the same in the purely geometrical case.

### Dynamic scintillation

Since scintillation is sensitive to small positional differences, it is also sensitive to small shifts of those over time. Because of the high velocity of most pulsars, scatterers can move quickly such that the source needs to be densely sampled within a short time to use the dynamic information. This can only be done with bright pulsars and rapidly repeating FRBs.

Velocities need to be assigned to the observer, the pulsar, as well as each scatterer. Velocity components along the line of sight are negligible because of the astronomical distances involved, such that all velocities are two-dimensional vectors:
~~~
D_s = 2100.*pc #distance to your favourite pulsar/FRB
V_p = [30.*kms,0.] #Earth's orbital velocity is 30 km/s
PMRA = -73.64
PMDEC = -62.65
V_RA = PMRA*mas/year*D_s
V_DEC = PMDEC*mas/year*D_s
V_s = [V_RA,V_DEC]
Sc = ScatteredSignal(D = 2100.*pc,V_p=V_p,V_s=V_s)
Sc.addPoint(D=125.*pc,x=0.1*au,V_x=1.*kms,V_y=0.1*kms)
~~~
After defining a time and a frequency axis, we can compute the dynamic spectrum, which is the same as above but with the movements of all objects taken into account:
~~~
t = np.linspace(0.,30.*minute,num=200)
nu = np.linspace(1200.*MHz,1600.*MHz,num=200)
I = Sc.compute_dynamic_spectrum(t,nu)
~~~
Two-dimensional data can be visualized with a colormap of your choice:
~~~
figure = plt.figure(figsize=(16,9))
ax = figure.add_subplot(1,1,1)
plot = ax.pcolormesh(t/minute,nu/MHz,np.swapaxes(I/np.std(I),0,1),cmap="viridis",vmin=0.,vmax=None,shading='nearest')
figure.colorbar(plot, ax=ax)
ax.set_xlabel(r"Time $t$ [minutes]")
ax.set_ylabel(r"frequency $\nu$ [MHz]")
~~~

~~~
Sec = np.abs(np.fft.fft2(I))**2
Sec = np.fft.fftshift(Sec)

N_nu = len(nu)
dnu = nu[1]-nu[0]
f_nu = np.fft.fftfreq(N_nu,dnu)
f_nu = np.fft.fftshift(f_nu)

N_t = len(t)
dt = t[1]-t[0]
f_t = np.fft.fftfreq(N_t,dt)
f_t = np.fft.fftshift(f_t)
~~~

~~~
figure = plt.figure(figsize=(16,9))
ax = figure.add_subplot(1,1,1)
plot = ax.pcolormesh(f_t/MHz,f_nu/mus,np.swapaxes(Sec,0,1),cmap="viridis",vmin=None,vmax=None,shading='nearest')
figure.colorbar(plot, ax=ax)
ax.set_xlabel(r"Doppler rate $f_{t}=f_{\rm D}$ [mHz]")
ax.set_ylabel(r"Delay $f_{\nu}=\tau$ [µs]")
~~~

~~~
Sec[Sec < 0.] = 0.
min_nonzero = np.min(Sec[np.nonzero(Sec)])
Sec[Sec == 0.] = min_nonzero
Sec = np.log10(Sec)
~~~

~~~
Sc.plot_delay_doppler(np.mean(nu))
~~~