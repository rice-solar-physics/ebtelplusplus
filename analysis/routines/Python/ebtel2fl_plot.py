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
        if 'parent_dir' not in kwargs:
            print "Warning: Parent directory not specified. You will not be able to load variables into class instance."
        else:
            self.parent_dir = kwargs['parent_dir']
        if 'alpha' in kwargs:
            self.alpha = kwargs['alpha']
        if 'tpulse' in kwargs:
            self.tpulse = kwargs['tpulse']
        if 'solver' in kwargs:
            self.solver = kwargs['solver']
        if 'looplength' in kwargs:
            self.looplength = kwargs['looplength']
        if 'species' in kwargs:
            self.species = kwargs['species']
        if 'Tn' in kwargs and 'run' in kwargs:
            self.Tn,self.run = kwargs['Tn'],kwargs['run']
            self.load_variables()
        elif 'data_file' in kwargs:
            self.load_variables(data_file=kwargs['data_file'])
        #configure static parameters
        self.fs = 18.0
        self.figsize = (10,10)
        self.linestyles = ('-','--','-.',':')
        self.format = 'eps'
        self.dpi = '1000'
            
            
    def load_variables(self,**kwargs):
        if 'data_file' in kwargs:
            data = np.loadtxt(self.parent_dir+kwargs['data_file'])
        else:
            child_dirs = self.species+'_heating_runs/alpha'+self.alpha+'/data/'
            child_file = 'ebtel2fl_L'+str(self.looplength)+'_tn'+str(self.Tn)+'_tpulse'+str(self.tpulse)+'_'+self.solver+'_'+str(self.run)+'.txt'
            data = np.loadtxt(self.parent_dir+child_dirs+child_file)
            
        self.time = data[:,0]
        self.temp_e = data[:,1]
        self.temp_i = data[:,2]
        self.dens = data[:,3]
        self.temp_apex_e = data[:,7]
        self.temp_apex_i = data[:,8]
        self.dens_apex = data[:,9]
        self.heat = data[:,15]
        
        
    def plot_em_curves(self,temp_list,em_list,**kwargs):
        #spacing between tn curves (artificial)
        delta_em = 0.2
        
        #set up figure
        fig = plt.figure(figsize=self.figsize)
        ax = fig.gca()
        
        #print lines
        for i in range(len(em_list)):
            if len(np.shape(np.array(em_list[i]))) > 1:
                mean_em = np.mean(em_list[i],axis=0)
                std_em = np.std(em_list[i],axis=0)
                mean_temp = np.mean(temp_list[i],axis=0)
                ax.plot(mean_temp,mean_em+i*delta_em,linestyle=self.linestyles[i%len(self.linestyles)],color='blue')
            else:
                ax.plot(temp_list[i],em_list[i]+i*delta_em,linestyle=self.linestyles[i%len(self.linestyles)],color='blue')
        
        #set labels
        ax.set_title(r'EBTEL Two-fluid EM, $\alpha$ = '+str(self.alpha),fontsize=self.fs)
        ax.set_xlabel(r'$\log T$ (K)',fontsize=self.fs)
        ax.set_ylabel(r'$\log$EM (cm$^{-5}$)',fontsize=self.fs)
        ax.set_xlim([5.5,7.5])
        ax.set_ylim([27,33])
        
        #save or show the figure
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename']+'.'+self.format,format=self.format,dpi=self.dpi)
        else:
            plt.show()
                