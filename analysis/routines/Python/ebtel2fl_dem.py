#ebtel2fl_dem.py

#Will Barnes
#13 May 2015

#Import needed modules
import numpy as np
import matplotlib.pyplot as plt

class DEMAnalyzer(object):
    
    def __init__(self,root_dir,species,alpha,loop_length,tpulse,solver,**kwargs):
        self.root_dir = root_dir
        self.species = species
        self.alpha = alpha
        self.loop_length = loop_length
        self.tpulse = tpulse
        self.solver = solver
        child_path = self.root_dir+self.species+'_heating_runs/alpha'+str(self.alpha)+'/data/'
        self.file_path = 'ebtel2fl_L'+str(self.loop_length)+'_tn%d_tpulse'+str(self.tpulse)+'_'+self.solver
        self.root_path = child_path + self.file_path
        if 'Tna' in kwargs:
            self.Tna = kwargs['Tna']
        else:
            self.Tna = 250
        if 'Tnb' in kwargs:
            self.Tnb = kwargs['Tnb']
        else:
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
            
            
    def process_raw(self,**kwargs):
        self.em = []
        self.temp_em = []
        for i in range(len(self.Tn)):
            tn_path = self.root_path%self.Tn[i]
            if self.mc is False:
                temp = np.loadtxt(tn_path+'_dem.txt')
                self.temp_em.append(temp[:,0])
                self.em.append(temp[:,4])
            else:
                em = []
                temp_em = []
                for j in range(self.mc):
                    try:
                        temp = np.loadtxt(tn_path+'/'+self.file_path%self.Tn[i]+'_'+str(j)+'_dem.txt')
                        temp_em.append(temp[:,0])
                        em.append(temp[:,4])
                    except:
                        raw_input("Unable to process file for Tn = "+str(self.Tn[i])+", run = "+str(j))
                        pass
                self.temp_em.append(temp_em)
                self.em.append(em)
                    