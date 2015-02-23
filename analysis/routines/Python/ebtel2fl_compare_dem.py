#Name: ebtel2fl_compare_dem.py
#Author: Will Barnes
#Date: 21 February 2014

#Import needed modules to plot non-interactively
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import argparse

#Declare parser object
parser = argparse.ArgumentParser(description='Script that collects DEM data from results of 2-fluid EBTEL model and shows coronal DEM and T_max plots.')

#Add arguments to parser
parser.add_argument("-s","--species",help="Species to which the heating was applied for particular run.")
parser.add_argument("-a","--alpha",type=float,help="Spectral index for the power-law distribution used.")
parser.add_argument("-L","--loop_length",type=float,help="Loop half-length.")
parser.add_argument("-t","--t_pulse",type=float,help="Width of the heating pulse used for the particular run.")

#Declare the parser dictionary
args = parser.parse_args()

#Set root directory for reading
species = args.species
root_dir = '/data/datadrive2/EBTEL-2fluid_runs/' + species + '_heating_runs/'
#Set the particular value of alpha
alpha = args.alpha
alpha_dir = 'alpha' + str(alpha) + '/' 
data_dir = root_dir + alpha_dir + 'data/'

#Set the loop length and pulse time
L = args.loop_length
t_pulse = args.t_pulse

#Create the array of wait times
wait_times = np.arange(250,5250,250)

#Set up the figure
fig1 = plt.figure(figsize=(10,10))
fig2 = plt.figure(figsize=(10,10))
ax1 = fig1.gca()
ax2 = fig2.gca()
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
    #Find the max value 
    ind_max = np.argmax(dem_cor)
    #Find the temperature at which the max occurs
    temp_max = tdem(ind_max)
    #Plot the values
    ax1.plot(tdem,dem_cor,linestyle=line_styles[i%4],color='blue')
    ax2.plot(wait_times(i),temp_max,'ko')

#Set some figure properties for the DEM plots
ax1.set_title(r'EBTEL Two-fluid DEM, $T_N=250-5000$ s',fontsize=fs)
ax1.set_xlabel(r'$\log(T_{DEM})$ (K)',fontsize=fs)
ax1.set_ylabel(r'$\log($DEM$)$ (cm$^{-5}$ K$^{-1}$)',fontsize=fs)
ax1.text(4.6,27.0,r'$\alpha$ = '+str(alpha),fontsize=fs)
ax1.set_xlim([4.5,7.0])
ax1.set_ylim([21,27.5])
#Set some properties for the Tmax plots
ax2.set_title(r'EBTEL Two-fluid $T(\max(DEM_C))$, $T_N=250-5000$ s',fontsize=fs)
ax2.set_xlabel(r'$T_N$',fontsize=fs)
ax2.set_ylabel(r'$\log(T_{max})$',fontsize=fs)
ax1.text(260,6.5,r'$\alpha$ = '+str(alpha),fontsize=fs)


#Save the figure to the top level directory
fig1
plt.savefig(root_dir+alpha_dir+'ebtel2fl_L'+str(L)+'_tpulse'+str(t_pulse)+'_alpha'+str(alpha)+ '_' + species + '_heating_dem.eps',format='eps',dpi=1000)
fig2
plt.savefig(root_dir+alpha_dir+'ebtel2fl_L'+str(L)+'_tpulse'+str(t_pulse)+'_alpha'+str(alpha)+ '_' + species + '_heating_TmaxVTn.eps',format='eps',dpi=1000)
    
