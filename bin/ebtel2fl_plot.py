#ebtel2fl_plot.py

#Will Barnes
#7 May 2015

#Import necessary modules
import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as aconst
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit

class Plotter(object):
    
    def __init__(self,**kwargs):
        #configure keyword arguments
        if 'parent_dir' in kwargs:
            self.parent_dir = parent_dir
        if 'child' in kwargs:
            self.child = child
        #configure static parameters
        self.fs = 18.0
        self.figsize = (10,10)
        self.linestyles = ('-','--','-.',':')
        self.format = 'eps'
        self.dpi = '1000'
        #load variables
        if 'parent_dir' in kwargs and 'child' in kwargs:
            self.load_variables()
        else:
            print "No file specified. Variable namespace will not be populated."    
            
    def load_variables(self,**kwargs):
        #load plasma parameters
        data = np.loadtxt(self.parent_dir+self.child+'.txt')
            
        self.time = data[:,0]
        self.temp_e = data[:,1]
        self.temp_i = data[:,2]
        self.dens = data[:,3]
        self.temp_apex_e = data[:,7]
        self.temp_apex_i = data[:,8]
        self.dens_apex = data[:,9]
        self.heat = data[:,15]
        
        #load dem parameters
        try:
            data = np.loadtxt(self.parent_dir+self.child+'_dem.txt')
            self.temp_dem = data[:,0]
            self.dem_tr = data[:,1]
            self.dem_cor = data[:,2]
            self.dem_tot = data[:,3]
            self.em_cor = data[:,4]
        except:
            print "Unable to load DEM parameters."
            pass
        #load heat parameters
        try:
            self.events = np.loadtxt(self.parent_dir+self.child+'_heat_amp.txt')
        except:
            print "Unable to load heating event amplitudes."
            pass
            
            
    def plot_params(self,**kwargs):
        #set up figure
        fig,ax = plt.subplots(3,1,figsize=(1.5*self.figsize[0],self.figsize[1]))
        ax_n = ax[1].twinx()
        ax_na = ax[2].twinx()
        
        #plot heating
        ax[0].plot(self.time,self.heat)
        ax[0].set_ylabel(r'$h$ (erg cm$^{-3}$ s$^{-1}$)',fontsize=self.fs)
        ax[0].set_title(r'EBTEL Two-fluid Plasma Parameters',fontsize=self.fs)
        ax[0].set_xlim([self.time[0],self.time[-1]])
        ax[0].locator_params(nbins=5)
        ax[0].ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
        #plot average temperature and density
        line_te = ax[1].plot(self.time,self.temp_e/10**6,label=r'$T_e$')
        line_ti = ax[1].plot(self.time,self.temp_i/10**6,'r--',label=r'$T_i$')
        ax[1].set_ylabel(r'$T$ (MK)',fontsize=self.fs)
        ax[1].yaxis.set_major_locator(MaxNLocator(prune='lower'))
        ax[1].locator_params(nbins=5)
        ax[1].ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
        line_n = ax_n.plot(self.time,self.dens/10**8,'k',label=r'$n$')
        ax_n.set_ylabel(r'$n$ (10$^8$ cm$^{-3}$)',fontsize=self.fs)
        ax_n.yaxis.set_major_locator(MaxNLocator(prune='lower'))
        ax_n.locator_params(nbins=5)
        ax_n.ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
        ax[1].set_xlim([self.time[0],self.time[-1]])
        #plot apex temperature and density
        ax[2].plot(self.time,self.temp_apex_e/10**6)
        ax[2].plot(self.time,self.temp_apex_i/10**6,'r--')
        ax[2].set_ylabel(r'$T_a$ (MK)',fontsize=self.fs)
        ax[2].yaxis.set_major_locator(MaxNLocator(prune='lower'))
        ax[2].locator_params(nbins=5)
        ax[2].ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
        ax_na.plot(self.time,self.dens_apex/10**8,'k')
        ax_na.set_ylabel(r'$n_a$ (10$^8$ cm$^{-3}$)',fontsize=self.fs)
        ax_na.yaxis.set_major_locator(MaxNLocator(prune='lower'))
        ax_na.locator_params(nbins=5)
        ax_na.ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
        ax[2].set_xlim([self.time[0],self.time[-1]])
        ax[2].set_xlabel(r'$t$ (s)',fontsize=self.fs)
        
        #configure legend
        lines = line_te + line_ti + line_n
        labels = [l.get_label() for l in lines]
        ax[1].legend(lines,labels,loc=1)
        
        #Check if output filename is specified
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename'],format=self.format,dpi=self.dpi)
        else:
            plt.show()
            
            
    def plot_dem(self,**kwargs):
        #set up figure
        fig = plt.figure(figsize=self.figsize)
        ax = fig.gca()
        
        #plot dem curves
        ax.plot(self.temp,self.dem_tr,label=r'TR')
        ax.plot(self.temp,self.dem_cor,label=r'corona')
        ax.plot(self.temp,self.dem_tot,label=r'total')
        ax.legend()
        ax.set_title(r'EBTEL Two-fluid DEM',fontsize=self.fs)
        ax.set_xlabel(r'$\log(T_{DEM})$ (K)',fontsize=self.fs)
        ax.set_ylabel(r'$\log($DEM$)$ (cm$^{-5}$ K$^{-1}$)',fontsize=self.fs)
        ax.set_xlim([5.5,7.5])

        #Check if output filename is specified
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename'],format=self.format,dpi=self.dpi)
        else:
            plt.show()
            
            
    def plot_event_distribution(self,**kwargs):
        #set up figure
        fig = plt.figure(figsize=self.figsize)
        ax = fig.gca()
        
        def power_law_curve(x,a,b):
            return a*(x**b)
        
        #Create a histogram and calculate fit
        dist,bins = np.histogram(self.events,bins=30)        
        pars,covar = curve_fit(power_law_curve,bins[0:-1],dist)
        pl_fit = power_law_curve(bins[0:-1],*pars)
        sigma = np.sqrt(np.diag(covar))

        #plot histogram
        ax.plot(bins[0:-1],dist,'ko',label=r'Events')
        ax.plot(bins[0:-1],pl_fit,'--r',label=r'Fit')
        ax.set_xlabel(r'Event Amplitude (erg cm$^{-3}$ s$^{-1}$)',fontsize=self.fs)
        ax.set_ylabel(r'Number of Events',fontsize=self.fs)
        ax.set_title(r'$P(x)=Cx^{\alpha}$, C = %.2e, $\alpha$ = %.2f $\pm$ %.2e' % (pars[0],pars[1],sigma[0]),fontsize=self.fs)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.legend(loc=1)
    
        #Check if output filename is specified
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename'],format=self.format,dpi=self.dpi)
        else:
            plt.show()
            
            
    def plot_surface(self,param_1,param_2,surf_list,**kwargs):
        #set up figure
        fig = plt.figure(figsize=self.figsize)
        ax = fig.gca()
        
        #set up mesh
        p1_mesh,p2_mesh = np.meshgrid(param_1,param_2)
        surf = ax.pcolor(p1_mesh,p2_mesh,np.array(surf_list),cmap='hot',vmin=np.min(np.array(surf_list)),vmax=np.max(np.array(surf_list)))
        fig.colorbar(surf,ax=ax)
        
        #set labels
        if 'ylab' in kwargs:
            ax.set_ylabel(kwargs['ylab'],fontsize=self.fs)
        if 'xlab' in kwargs:
            ax.set_xlabel(kwargs['xlab'],fontsize=self.fs)
        if 'plot_title' in kwargs:
            ax.set_title(kwargs['plot_title'],fontsize=self.fs)
            
        #Check if output filename is specified
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename'],format=self.format,dpi=self.dpi)
        else:
            plt.show()
