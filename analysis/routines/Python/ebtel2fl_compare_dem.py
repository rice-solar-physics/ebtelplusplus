#Name: ebtel2fl_compare_dem.py
#Author: Will Barnes
#Date: 21 February 2014

#Import needed modules to plot non-interactively
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

#Set root directory for reading
species = 'ion'
root_dir = '/data/datadrive2/EBTEL-2fluid_runs/' + species + '_heating_runs/'
#Set the particular value of alpha
alpha = 1.5
alpha_dir = 'alpha' + str(alpha) + '/' 
data_dir = root_dir + alpha_dir + 'data/'

#Set the loop length and pulse time
L = 20.0
t_pulse = 200.0

#Create the array of wait times
wait_times = np.arange(250,5250,250)

#Set up the figure
#Set up the figure
fig = plt.figure(figsize=(10,10))
ax = fig.gca()
fs = 18

#Set linestyle options
line_styles = ('-','--','-.',':')
#Set artificial spacing between lines
delta = 0.3

#Start the loop to read in the values
for i in range(len(wait_times)):
    #Make the dir name
    temp_file = data_dir + 'ebtel2fl_L' + str(L) + '_tn' + str(wait_times[i]) + '_tpulse' + str(t_pulse) + '_dem.txt'
    #Load the text file
    temp = np.loadtxt(temp_file)
    #Get the logTdem and dem_cor values
    tdem = temp[:,0]
    dem_cor = temp[:,2] + i*delta
    #Plot the values
    ax.plot(tdem,dem_cor,linestyle=line_styles[i%4],color='blue')
    

#Set some figure properties
ax.set_title(r'EBTEL Two-fluid DEM, $T_N=250-5000$ s',fontsize=fs)
ax.set_xlabel(r'$\log(T_{DEM})$ (K)',fontsize=fs)
ax.set_ylabel(r'$\log($DEM$)$ (cm$^{-5}$ K$^{-1}$)',fontsize=fs)
ax.text(4.6,27.0,r'$\alpha$ = '+str(alpha),fontsize=fs)
ax.set_xlim([4.5,7.0])
ax.set_ylim([21,27.5])

#Save the figure to the top level directory
plt.savefig(root_dir+alpha_dir+'ebtel2fl_L'+str(L)+'_tpulse'+str(t_pulse)+'_alpha'+str(alpha)+ '_' + species + '_heating_dem.eps',format='eps',dpi=1000)
    
