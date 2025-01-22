import numpy as np
from numpy import newaxis as na
import matplotlib as mpl
import matplotlib.pyplot as plt

def draw_canvas(plot_width = 1200, plot_height = 900, plot_dpi = 100, plot_bottom = 0.08, plot_top = 0.95, plot_left = 0.08, plot_right = 0.98, plot_wspace = 0.1, plot_hspace = 0.2, textsize=20, labelsize=18):
    figure = plt.figure(figsize=(plot_width/plot_dpi,plot_height/plot_dpi),dpi=plot_dpi)
    plt.subplots_adjust(bottom=plot_bottom,top=plot_top,left=plot_left,right=plot_right,wspace=plot_wspace,hspace=plot_hspace)
    pgf_with_pdflatex = {
        "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
        "text.usetex": True,                # use LaTeX to write all text
        "font.size": labelsize,
        "axes.linewidth": 0.5,                
        "axes.labelsize": labelsize,               # LaTeX default is 10pt font. 
        "axes.titlesize": labelsize,
        "patch.linewidth": 0.5,		# Width of box around legend
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "lines.markeredgewidth": 0.3,
        "legend.fontsize": textsize,
        "legend.edgecolor": "black",
        "legend.borderpad": 0.3,			# width of whitespace between text and border (in units of times linewidth)
        "xtick.labelsize": textsize,
        "ytick.labelsize": textsize,
    }
    mpl.rcParams.update(pgf_with_pdflatex) 
    
    return figure

#constants
au = 149597870700. #m
pc = 648000./np.pi*au #m
v_c = 299792458. #m/s
ms = 1.0e-3 #s
mus = 1.0e-6 #s
degrees = np.pi/180.
mas = degrees/1000./3600.
MHz = 1.0e+6 #Hz
kms = 1000. #km/s

class ScatteredSignal:
    def __init__(self,**kwargs):
        self.D_s = kwargs.get("D",1000.*pc)
        self.V_s = np.array(kwargs.get("V_s",[0.,0.]))
        self.V_p = np.array(kwargs.get("V_p",[0.,0.]))
        
        #initialize scattered points with central los
        self.points = np.empty((6,0),dtype=float)
        
    def addPoint(self,**kwargs):
        x = kwargs.get("x",1.*au)
        y = kwargs.get("y",0.)
        D_x = kwargs.get("D",self.D_s/2.)
        mu = kwargs.get("mu",0.5)
        V_x = kwargs.get("V_x",0.)
        V_y = kwargs.get("V_y",0.)
        point = np.array([x,y,D_x,mu,V_x,V_y])
        
        self.points = np.concatenate((self.points,point[:,na]),axis=1)
        
    def plot_rays(self,**kwargs):
        figure = draw_canvas(plot_width = 1200,plot_height = 700)
        ax = figure.add_subplot(1,1,1)
        ax.plot([0.,self.D_s/pc],[0.,0.],color='black',marker="",linestyle="-")
        ax.set_xlabel("Distance along line of sight [pc]")
        ax.set_ylabel("Distance perpendicular to los [au]")
        for i_point in range(self.points.shape[1]):
            ax.plot([0.,self.points[2,i_point]/pc,self.D_s/pc],[0.,self.points[0,i_point]/au,0.],marker="",linestyle="-",color="blue")
            
    def plot_delays(self,**kwargs):
        Deff = self.points[2]*self.D_s/(self.D_s-self.points[2])
        delays = Deff/(2.*v_c)*(self.points[0]**2+self.points[1]**2)/self.points[2]**2
        
        figure = draw_canvas(plot_width = 1200,plot_height = 700)
        ax = figure.add_subplot(1,1,1)
        ax.vlines([0.],[0.],[1.],color='black')
        ax.set_xlabel(r"Delay $\tau$ [µs]")
        ax.set_ylabel("Intensity $I$ (1=intrinsic)")
        ax.set_ylim([0.,None])
        for i_point in range(self.points.shape[1]):
            ax.vlines(delays[i_point]/mus,[0],self.points[3,i_point],color="blue")
        #ax.vlines(delays/ms,np.zeros(self.points.shape[1]),self.points[3],color="blue")
        
    def plot_images(self,**kwargs):
        thx = self.points[0]/self.points[2]
        thy = self.points[1]/self.points[2]
        
        figure = draw_canvas(plot_width = 1200,plot_height = 700)
        ax = figure.add_subplot(1,1,1)
        ax.set_xlabel(r"$\theta_x$ [mas]")
        ax.set_ylabel(r"$\theta_y$ [mas]")
        ax.plot([0],[0],color="black",marker='o',linestyle="")
        ax.plot(thx/mas,thy/mas,color="blue",marker='o',linestyle="")
        
    def compute_scintillation(self,nu,**kwargs):
        Deff = self.points[2]*self.D_s/(self.D_s-self.points[2])
        delays = Deff/(2.*v_c)*(self.points[0]**2+self.points[1]**2)/self.points[2]**2
        
        E = np.ones_like(nu,dtype=complex)
        for i_nu,v_nu in enumerate(nu):
            E[i_nu] += np.sum(self.points[3]*np.exp(2.0j*np.pi*v_nu*delays))
        return np.abs(E)**2
        
    def plot_scintillation(self,nu,I,**kwargs):
        figure = draw_canvas(plot_width = 1200,plot_height = 700)
        ax = figure.add_subplot(1,1,1)
        ax.set_xlabel(r"frequency $\nu$ [MHz]")
        ax.set_ylabel(r"Intensity $I$ (max=1)")
        ax.set_ylim([0.,1.1])
        ax.plot(nu/MHz,I/np.max(I),color="black",marker='',linestyle="-")
        
    def plot_differential_delays(self,**kwargs):
        y_max = kwargs.get("y_max",1.)
        
        Deff = self.points[2]*self.D_s/(self.D_s-self.points[2])
        delays = Deff/(2.*v_c)*(self.points[0]**2+self.points[1]**2)/self.points[2]**2
        
        figure = draw_canvas(plot_width = 1200,plot_height = 700)
        ax = figure.add_subplot(1,1,1)
        ax.set_xlabel(r"Delay $\tau$ [µs]")
        ax.set_ylabel("Intensity $I$ (1=center)")
        ax.set_ylim([0.,y_max])
        mu_center = (1.+np.sum(self.points[3]**2))**2
        for i_point in range(self.points.shape[1]):
            ax.vlines(delays[i_point]/mus,[0],self.points[3,i_point]**2/mu_center,color="blue")
            ax.vlines(-delays[i_point]/mus,[0],self.points[3,i_point]**2/mu_center,color="blue")
            for i_point2 in range(self.points.shape[1]):
                if i_point!=i_point2:
                    ax.vlines((delays[i_point]-delays[i_point2])/mus,[0],(self.points[3,i_point]**2*self.points[3,i_point2]**2)/mu_center,color="red")
        ax.vlines([0.],[0.],[1.],color='black')
        
    def compute_dynamic_spectrum(self,t,nu,**kwargs):
        Deff = self.points[2]*self.D_s/(self.D_s-self.points[2])
        delays = Deff/(2.*v_c)*(self.points[0]**2+self.points[1]**2)/self.points[2]**2
        
        Veff = self.V_p[:,na] - self.D_s/(self.D_s-self.points[2])*self.points[(4,5),:] + self.points[2]/(self.D_s-self.points[2])*self.V_s[:,na]
        dopplers = -1.*(Veff[0]*self.points[0]+Veff[1]*self.points[1])/(self.points[2]*v_c)
        
        E = np.ones((len(t),len(nu)),dtype=complex)
        for i_nu,v_nu in enumerate(nu):
            for i_t,v_t in enumerate(t):
                E[i_t,i_nu] += np.sum(self.points[3]*np.exp(2.0j*np.pi*v_nu*(delays+v_t*dopplers)))
        return np.abs(E)**2