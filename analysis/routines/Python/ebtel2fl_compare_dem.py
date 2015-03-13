#Name: ebtel2fl_compare_dem.py
#Author: Will Barnes
#Date: 21 February 2014

#Import needed modules to plot non-interactively
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_ebtel_dem_compare(species,alpha,L,t_pulse,solver,**kwargs):
    """Plot the DEM for different runs of EBTEL-2fluid with differing waiting times.

    Arguments:
    species -- which species is being heated
    alpha -- spectral power-law index (>0 here)
    t_pulse -- duration of heating pulse (s)
    solver -- solver option being used

    Optional keyword arguments:
    print_fig_filename -- set to 'true' to print to file; recommended for large number of plots.

    """
    
    #Set root directory for reading
    root_dir = '/data/datadrive2/EBTEL-2fluid_runs/' + species + '_heating_runs/'
    alpha_dir = 'alpha' + str(alpha) + '/'
    data_dir = root_dir + alpha_dir + 'data/'

    #Create the array of wait times
    wait_times = np.arange(250,5250,250)

    #Set up the figure
    fig1 = plt.figure(figsize=(10,10))
    fig2 = plt.figure(figsize=(10,7))
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
        temp_file = data_dir + 'ebtel2fl_L' + str(L) + '_tn' + str(wait_times[i]) + '_tpulse' + str(t_pulse) + '_'+ solver + '_dem.txt'
        #Load the text file
        temp = np.loadtxt(temp_file)
        #Get the logTdem and dem_cor values
        tdem = temp[:,0]
        dem_cor = temp[:,2] + i*delta
        #Find the max value
        ind_max = np.argmax(dem_cor)
        #Find the temperature at which the max occurs
        temp_max = tdem[ind_max]
        #Plot the values
        ax1.plot(tdem,dem_cor,linestyle=line_styles[i%4],color='blue')
        ax2.plot(wait_times[i],temp_max,'ko')

    #Set some figure properties for the DEM plots
    ax1.set_title(r'EBTEL Two-fluid DEM, $T_N=250-5000$ s',fontsize=fs)
    ax1.set_xlabel(r'$\log(T_{DEM})$ (K)',fontsize=fs)
    ax1.set_ylabel(r'$\log($DEM$)$ (cm$^{-5}$ K$^{-1}$)',fontsize=fs)
    ax1.text(4.6,27.0,r'$\alpha$ = '+str(alpha),fontsize=fs)
    ax1.set_xlim([4.5,7.5])
    ax1.set_ylim([22,28.5])
    #Set some properties for the Tmax plots
    ax2.set_title(r'EBTEL Two-fluid $T(\max(DEM_C))$, $T_N=250-5000$ s',fontsize=fs)
    ax2.set_xlabel(r'$T_N$',fontsize=fs)
    ax2.set_ylabel(r'$\log(T_{max})$',fontsize=fs)
    ax2.text(500,6.8,r'$\alpha$ = '+str(alpha),fontsize=fs)
    ax2.set_ylim([5.5,7.0])
    
    #Check if output filename is specified
    if 'print_fig_filename' in kwargs:
        #Save the figures
        plt.figure(fig1.number)
        plt.savefig(root_dir+alpha_dir+'ebtel2fl_L'+str(L)+'_tpulse'+str(t_pulse)+'_alpha'+str(alpha)+ '_' + species + '_heating_dem.eps',format='eps',dpi=1000)
        plt.figure(fig2.number)
        plt.savefig(root_dir+alpha_dir+'ebtel2fl_L'+str(L)+'_tpulse'+str(t_pulse)+'_alpha'+str(alpha)+ '_' + species + '_heating_TmaxVTn.eps',format='eps',dpi=1000)
    else:
        plt.figure(fig1.number)
        plt.show()
        plt.figure(fig2.number)
        plt.show()


