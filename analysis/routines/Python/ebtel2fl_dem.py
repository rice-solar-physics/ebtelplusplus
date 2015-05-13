#ebtel2fl_dem.py

#Will Barnes
#13 May 2015

#Import needed modules
import numpy as np
import matplotlib.pyplot as plt

class DEMAnalyzer(object):
    
    def __init__(self,alpha,loop_length,tpulse,solver,**kwargs):
        self.alpha = alpha
        self.loop_length = loop_length
        self.tpulse = tpulse
        self.solver = solver
        if 'Tna' in kwargs:
            self.Tna = kwargs['Tna']
        else:
            self.Tna = 250
        if 'Tnb' in kwargs:
            self.Tnb = kwargs['Tnb']
        else
            self.Tnb = 5000
        if 'Tndelta' in kwargs:
            self.Tndelta = kwargs['Tndelta']
        else:
            self.Tndelta = 250
        self.Tn = np.arange(self.Tna,self.Tnb+self.Tndelta,self.Tndelta)
        if 'mc' in kwargs:
            self.mc = kwargs['mc']
        else:
            self.mc = False