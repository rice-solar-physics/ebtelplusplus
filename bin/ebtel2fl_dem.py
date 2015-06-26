#ebtel2fl_dem.py

#Will Barnes
#13 May 2015

#Import needed modules
import numpy as np
from scipy.optimize import curve_fit


class DEMAnalyzer(object):
    
    def __init__(self,root_dir,species,alpha,loop_length,tpulse,solver,**kwargs):
        #set object variables
        self.root_dir = root_dir
        self.species = species
        self.alpha = alpha
        self.loop_length = loop_length
        self.tpulse = tpulse
        self.solver = solver
        #set up paths
        child_path = self.root_dir+self.species+'_heating_runs/alpha'+str(self.alpha)+'/data/'
        self.file_path = 'ebtel_L'+str(self.loop_length)+'_tn%d_tpulse'+str(self.tpulse)+'_'+self.solver
        self.root_path = child_path + self.file_path
        #configure keyword arguments
        if 'Tn' in kwargs:
            self.Tn = kwargs['Tn']
        else:
            self.Tn = np.arange(250,5250,250)
        if 'delta_cool' in kwargs:
            self.delta_cool = kwargs['delta_cool']
        else:
            self.delta_cool = -0.6
        if 'delta_hot' in kwargs:
            self.delta_hot = kwargs['delta_hot']
        else:
            self.delta_hot = -2.0
        #set static variables
        self.em_cutoff = 23.0
        self.em_max_eps_percent = 0.999
            
            
    def process_raw(self,**kwargs):
        self.em = []
        self.temp_em = []
        for i in range(len(self.Tn)):
            tn_path = self.root_path%self.Tn[i]
            #initialize lists
            em = []
            temp_em = []
            #initialize flag and counter
            eol_flag=False
            counter=0
            while eol_flag is False:
                try:
                    #load data
                    temp = np.loadtxt(tn_path+'/'+self.file_path%self.Tn[i]+'_'+str(counter)+'_dem.txt')
                    #check for nan
                    temp[np.where(np.isnan(temp))] = -np.inf
                    #append temperature and EM
                    temp_em.append(temp[:,0])
                    em.append(temp[:,4])
                    #increment counter
                    counter += 1
                except:
                    print "Unable to process file for Tn = "+str(self.Tn[i])+", run = "+str(counter)
                    print "Reached end of list or there was an error reading the file."
                    eol_flag=True
                    pass
            self.temp_em.append(temp_em)
            self.em.append(em)
                
                
    def em_max(self,**kwargs):
        self.temp_max = []
        self.em_max = []
        for i in range(len(self.Tn)):
            temp_max_temp = []
            em_max_temp = []
            for j in range(len(self.em[i])):
                i_max = np.argmax(self.em[i][j])
                indices_em_max = np.where(np.array(self.em[i][j]) > self.em_max_eps_percent*self.em[i][j][i_max])
                if indices_em_max <= 1:
                    temp_max_temp.append(self.temp_em[i][j][i_max])
                    em_max_temp.append(self.em[i][j][i_max])
                else:
                    em_interp = np.linspace(self.em[i][j][indices_em_max[0][0]],self.em[i][j][indices_em_max[0][-1]],100)
                    temp_interp = np.interp(em_interp,np.array(self.em[i][j][indices_em_max[0][0]:indices_em_max[0][-1]]),np.array(self.temp_em[i][j][indices_em_max[0][0]:indices_em_max[0][-1]]))
                    temp_max_temp.append(np.mean(temp_interp))
                    em_max_temp.append(np.mean(em_interp))
            self.temp_max.append(temp_max_temp)
            self.em_max.append(em_max_temp)
                
                
    def many_slopes(self,**kwargs):
        self.a_cool = []
        self.a_hot = []
        for i in range(len(self.Tn)):
            acl = []
            ahl = []
            for j in range(len(self.temp_em[i])):
                ac,ah = self.slope(self.temp_em[i][j],self.em[i][j])
                acl.append(ac),ahl.append(ah)
            self.a_cool.append(acl),self.a_hot.append(ahl)


    def slope(self,temp,dem,**kwargs):
        #Calculate bounds
        dict_bounds = self.bounds(temp,dem)
        
        #Function for linear fit
        def linear_fit(x,a,b):
            return a*x + b
            
        #Check if inside interpolated array and calculate slope
        #cool
        if np.size(dict_bounds['bound_cool']) == 0:
            print "Cool bound out of range. T_cool_bound = ",temp[np.argmax(dem)] + self.delta_cool," < T_cool(0) = ",dict_bounds['temp_cool'][0]
            a_coolward = False
        else:
            bound_cool = dict_bounds['bound_cool'][0][-1] + 1
            pars_cool,covar = curve_fit(linear_fit,dict_bounds['temp_cool'][bound_cool:-1],dict_bounds['dem_cool'][bound_cool:-1])
            a_coolward = pars_cool[0]
        #hot
        if np.size(dict_bounds['bound_hot']) == 0:
            print "Hot bound out of range. DEM_hot_bound = ",dem[np.argmax(dem)] + self.delta_hot," < DEM_hot(end) = ",dict_bounds['dem_hot'][-1]
            a_hotward = False
        else:
            bound_hot = dict_bounds['bound_hot'][0][0] - 1
            pars_hot,covar = curve_fit(linear_fit,dict_bounds['temp_hot'][0:bound_hot],dict_bounds['dem_hot'][0:bound_hot])
            a_hotward = pars_hot[0]
            
        return a_coolward,a_hotward
        
        
    def integrate(self,temp,dem,**kwargs):
        #Find the corresponding temperature bounds
        dict_bounds = self.bounds(temp,dem)

        #First check if the bounds are inside of our interpolated array
        if np.size(dict_bounds['bound_cool']) == 0 or np.size(dict_bounds['bound_hot']) == 0:
            print "Cool and/or hot bound(s) out of range. Skipping integration for these bounds."
            hot_shoulder_strength = False
        else:
            #Refine the arrays we will integrate over
            #Temprature
            temp_hot = dict_bounds['temp_hot'][0:(dict_bounds['bound_hot'][0][0] - 1)]
            temp_cool = dict_bounds['temp_cool'][(dict_bounds['bound_cool'][0][-1] + 1):-1]
            #DEM (EM)
            dem_hot = dict_bounds['dem_hot'][0:(dict_bounds['bound_hot'][0][0] - 1)]
            dem_cool = dict_bounds['dem_cool'][(dict_bounds['bound_cool'][0][-1] + 1):-1]
            #Do the integration
            #Hot shoulder
            hot_shoulder = np.trapz(dem_hot,x=temp_hot)
            #Total
            total_shoulder = np.trapz(np.concatenate([dem_cool[0:-1],dem_hot]),x=np.concatenate([temp_cool[0:-1],temp_hot]))
            #Compute the ratio
            hot_shoulder_strength = hot_shoulder/total_shoulder

        return hot_shoulder_strength
        
        
    def bounds(self,temp,dem,**kwargs):
        #Find peak parameters
        dem_max = np.max(dem)
        i_dem_max = np.argmax(dem)
        temp_dem_max = temp[i_dem_max]

        #Create cool and hot DEM and temperature arrays
        dem_hot = dem[i_dem_max:-1]
        temp_hot = temp[i_dem_max:-1]
        dem_cool = dem[0:i_dem_max]
        temp_cool = temp[0:i_dem_max]

        #Find the dem index where dem->inf (or less than the cutoff)
        inf_index_hot = np.where(dem_hot > self.em_cutoff)[0][-1]
        inf_index_cool = np.where(dem_cool > self.em_cutoff)[0][0]

        #Calculate the cool and hot bounds (in DEM and temperature)
        temp_cool_bound = temp_dem_max + self.delta_cool
        dem_hot_bound = dem_max + self.delta_hot

        #Interpolate and find accurate bound index
        temp_cool_new = np.linspace(temp_cool[inf_index_cool],temp_cool[-1],1000)
        dem_cool_new = np.interp(temp_cool_new,temp_cool[inf_index_cool:-1],dem_cool[inf_index_cool:-1])
        i_bound_cool = np.where(temp_cool_new < temp_cool_bound)
        temp_hot_new = np.linspace(temp_hot[0],temp_hot[inf_index_hot],1000)
        dem_hot_new = np.interp(temp_hot_new,temp_hot[0:inf_index_hot],dem_hot[0:inf_index_hot])
        i_bound_hot = np.where(dem_hot_new < dem_hot_bound)

        #Return the list of indices and interpolated DEM and temperature arrays
        return {'bound_cool':i_bound_cool,'bound_hot':i_bound_hot,'temp_cool':temp_cool_new,'temp_hot':temp_hot_new,'dem_cool':dem_cool_new,'dem_hot':dem_hot_new}
        