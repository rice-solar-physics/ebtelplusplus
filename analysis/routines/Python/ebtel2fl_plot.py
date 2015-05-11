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
    
    def __init__(self,parent_dir,**kwargs):
        self.parent_dir = parent_dir
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
        if 'data_file' in kwargs and kwargs is False:
            self.load_variables(data_file=kwargs['data_file'])
            
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