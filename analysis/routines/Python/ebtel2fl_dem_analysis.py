#Name: ebtel2fl_compare_dem.py
#Author: Will Barnes
#Date: 21 February 2014

#Import needed modules to plot non-interactively
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def find_temp_bounds(temp,dem,delta_cool,delta_hot):
    """Calculate corresponding temperature bounds for DEM threshold value.

    Arguments:
    temp -- log of temperature bin
    dem -- log of coronal DEM value
    delta_cool -- orders of magnitude below the DEM peak to set the cool bound
    delta_hot -- orders of magnitude below the DEM peak to set the hot bound; may need to be lower than delta_cool for ion heating
    
    """

    #Find peak DEM value
    dem_max = np.max(dem)

    #Find temperature for peak DEM value
    i_dem_max = np.argmax(dem)
    temp_dem_max = temp[i_dem_max]

    #Create cool and hot DEM and temperature arrays
    dem_hot = dem[i_dem_max:-1]
    temp_hot = temp[i_dem_max:-1]
    dem_cool = dem[0:i_dem_max]
    temp_cool = temp[0:i_dem_max]

    #Find the dem index where dem->inf for the hot side
    inf_index_hot = np.where(np.isinf(dem_hot)==False)[0][-1]
    #Find the dem index where dem->inf for the cool side
    inf_index_cool = np.where(np.isinf(dem_cool)==False)[0][0]

    #Calculate the cool and hot bounds (in DEM and temperature)
    #Cool shoulder
    temp_cool_bound = temp_dem_max - delta_cool
    #Hot shoulder
    temp_hot_bound = temp_dem_max + delta_hot

    #Check if our bounds are valid for these temp and dem arrays
    #If they are valid, calculate the hotward and coolward slopes
    #Interpolate over the cool branch
    temp_cool_new = np.linspace(temp_cool[inf_index_cool],temp_cool[-1],1000)
    dem_cool_new = np.interp(temp_cool_new,temp_cool[inf_index_cool:-1],dem_cool[inf_index_cool:-1])
    #Find the more accurate index of the cool bound
    i_bound_cool = np.where(temp_cool_new < temp_cool_bound)
    
    #Interpolate over the hot branch
    temp_hot_new = np.linspace(temp_hot[0],temp_hot[inf_index_hot],1000)
    dem_hot_new = np.interp(temp_hot_new,temp_hot[0:inf_index_hot],dem_hot[0:inf_index_hot])
    #Find the more accurate index of the hot bound
    i_bound_hot = np.where(temp_hot_new > temp_hot_bound)
    
    #Return the list of indices and interpolated DEM and temperature arrays
    return {'bound_cool':i_bound_cool,'bound_hot':i_bound_hot,'temp_cool':temp_cool_new,'temp_hot':temp_hot_new,'dem_cool':dem_cool_new,'dem_hot':dem_hot_new}
    
    
def dem_shoulder_compare_fit(temp,dem,delta_cool,delta_hot):
    """Compute coolward and hotward slope of DEM curve for a linear fit.

    Arguments:
    temp -- log of temperature bin
    dem -- log of coronal DEM value
    delta_cool -- orders of magnitude below the DEM peak to set the cool bound
    delta_hot -- orders of magnitude below the DEM peak to set the hot bound; may need to be lower than delta_cool for ion heating

    """
    
    #Call function to interpolate and find appropriate bounds
    dict_bounds = find_temp_bounds(temp,dem,delta_cool,delta_hot)
    
    #Declare function for our curve fit
    def linear_fit(x,a,b):
        return a*x + b
    
    #Check if the bound is inside of our interpolated array
    if np.size(dict_bounds['bound_cool']) == 0:
        print "Cool bound out of range. T_cool_bound = ",temp[np.argmax(dem)] - delta_cool," < T_cool(0) = ",dict_bounds['temp_cool'][0]
        a_coolward = False
    else:
        bound_cool = dict_bounds['bound_cool'][0][-1] + 1
        #Calculate the coolward slope
        pars_cool,covar = curve_fit(linear_fit,dict_bounds['temp_cool'][bound_cool:-1],dict_bounds['dem_cool'][bound_cool:-1])
        a_coolward = pars_cool[0]

    #Check if the bound is inside of the interpolated array
    if np.size(dict_bounds['bound_hot']) == 0:
        print "Hot bound out of range. T_hot_bound = ",temp[np.argmax(dem)] + delta_hot," > T_hot(end) = ",dict_bounds['temp_hot'][-1]
        a_hotward = False
    else:
        bound_hot = dict_bounds['bound_hot'][0][0] - 1
        #Calculate the hotward slope
        pars_hot,covar = curve_fit(linear_fit,dict_bounds['temp_hot'][0:bound_hot],dict_bounds['dem_hot'][0:bound_hot])
        a_hotward = pars_hot[0]

    #Return the hot and cool slopes
    return {'a_hot':a_hotward,'a_cool':a_coolward}



def dem_shoulder_compare_integrate(temp,dem,delta_cool,delta_hot):
    """Compute integral of hot and cold shoulder and calculate ratio to provide quantitative measure of hot DEM component.

    Arguments:
    temp -- log of temperature bin
    dem -- log of coronal DEM value
    delta_cool -- orders of magnitude below the DEM peak to set the cool bound
    delta_hot -- orders of magnitude below the DEM peak to set the hot bound; may need to be lower than delta_cool for ion heating

    """
    #Find the corresponding temperature bounds
    dict_bounds = find_temp_bounds(temp,dem,delta_cool,delta_hot)
    
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


def plot_ebtel_dem_compare(species,alpha,L,t_pulse,solver):
    """Plot the DEM (or EM) for different runs of EBTEL-2fluid with differing waiting times.

    Arguments:
    species -- which species is being heated
    alpha -- spectral power-law index (>0 here)
    t_pulse -- duration of heating pulse (s)
    solver -- solver option being used

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
    fig3,ax3 = plt.subplots(3,1,figsize=(12,10))
    ax1 = fig1.gca()
    ax2 = fig2.gca()
    ax2_em = ax2.twinx()
    fs = 18

    #Set linestyle options
    line_styles = ('-','--','-.',':')
    #Set artificial spacing between lines
    delta = 0.2

    #Start the loop to read in the values
    for i in range(len(wait_times)):
        #Make the dir name
        temp_file = data_dir + 'ebtel2fl_L' + str(L) + '_tn' + str(wait_times[i]) + '_tpulse' + str(t_pulse) + '_'+ solver + '_dem.txt'
        #Load the text file
        temp = np.loadtxt(temp_file)
        #Get the logTdem and dem_cor values
        tdem = temp[:,0]
        dem_cor = temp[:,4]
        #Find the max value
        ind_max = np.argmax(dem_cor)
        #Find the temperature at which the max occurs
        temp_max = tdem[ind_max]
        #Calculate the hot shoulder value and the fits; take an average over several different bounds for the cool and hot shoulders
        intervals_cool = np.linspace(0.6,0.6,1)
        intervals_hot = np.linspace(0.4,0.4,1)
        #Initialize integration and fits
        a_cool = []
        a_hot = []
        hs_int = []
        for j in range(len(intervals)):
            #Do the shoulder fit
            em_fit=dem_shoulder_compare_fit(tdem,dem_cor,intervals_cool[j],intervals_hot[j])
            if em_fit['a_cool'] != False:
                a_cool.append(em_fit['a_cool'])
            if em_fit['a_hot'] != False:
                a_hot.append(em_fit['a_hot'])
            #Do the integration
            em_int = dem_shoulder_compare_integrate(tdem,dem_cor,intervals_cool[j],intervals_hot[j])
            if  em_int != False:
                hs_int.append(em_int)

        #Average all of the values
        a_cool = np.mean(a_cool)
        a_hot = np.mean(a_hot)
        hs_int = np.mean(hs_int)
        #Plot the DEM (EM) values, adding an arbitrary separation
        ax1.plot(tdem,dem_cor+ i*delta,linestyle=line_styles[i%4],color='blue')
        #Plot the Tmax values
        ax2.plot(wait_times[i],temp_max,'ko')
        #Plot the EMmax values
        ax2_em.plot(wait_times[i],dem_cor[ind_max],'k+')
        #Plot the different shoulder strength measurements
        #Hot shoulder integration
        if np.isnan(hs_int) == False:
            ax3[0].plot(wait_times[i],hs_int,'ko')
        #EM(T) slope fits
        if np.isnan(a_hot) == False:
            ax3[1].plot(wait_times[i],abs(a_hot),'ro')
        if np.isnan(a_cool) == False:
            ax3[1].plot(wait_times[i],abs(a_cool),'bo')
        if np.isnan(a_hot) == False and np.isnan(a_cool) == False:
            ax3[2].plot(wait_times[i],abs(a_cool/a_hot),'ko')

    #Set some figure properties for the DEM plots
    ax1.set_title(r'EBTEL Two-fluid EM, $\alpha$ = '+str(alpha),fontsize=fs)
    ax1.set_xlabel(r'$\log T$ (K)',fontsize=fs)
    ax1.set_ylabel(r'$\log$EM (cm$^{-5}$)',fontsize=fs)
    ax1.set_xlim([5.5,7.5])
    ax1.set_ylim([27,33])
    #Set some properties for the Tmax plots
    ax2.set_title(r'EBTEL Two-fluid $T(\max(EM))$, $\alpha$ = '+str(alpha),fontsize=fs)
    ax2.set_xlabel(r'$T_N$',fontsize=fs)
    ax2.set_ylabel(r'$\log(T_{max})$',fontsize=fs)
    ax2.set_ylim([5.5,7.0])
    ax2.set_xlim([wait_times[0]-250,wait_times[-1]+250])
    #Set some properties for the EM max plots
    ax2_em.set_ylabel(r'$\log$EM($T_{max}$) (cm$^{-5}$)',fontsize=fs)
    ax2_em.set_ylim([26,30])
    #Set some properties for the hot shoulder strength comparison plots
    ax3[0].set_title(r'EBTEL Two-fluid Hot Shoulder Strength Comparison',fontsize=fs)
    ax3[0].set_ylabel(r'Integration',fontsize=fs)
    ax3[0].plot([wait_times[0]-250,wait_times[-1]+250],[0.5,0.5],'--k')
    ax3[0].set_ylim([0,1])
    ax3[0].set_xlim([wait_times[0]-250,wait_times[-1]+250])
    ax3[1].set_ylabel(r'$a_{hot,cool}$',fontsize=fs)
    ax3[1].plot([wait_times[0]-250,wait_times[-1]+250],[2,2],'--k')
    ax3[1].plot([wait_times[0]-250,wait_times[-1]+250],[3,3],'-k')
    ax3[1].plot([wait_times[0]-250,wait_times[-1]+250],[5,5],'-.k')
    ax3[1].set_ylim([0,10])
    ax3[1].set_xlim([wait_times[0]-250,wait_times[-1]+250])
    ax3[2].set_ylabel(r'$a_{cool}/a_{hot}$',fontsize=fs)
    ax3[2].plot([wait_times[0]-250,wait_times[-1]+250],[1,1],'--k')
    ax3[2].set_ylim([0,5])
    ax3[2].set_xlim([wait_times[0]-250,wait_times[-1]+250])
    ax3[2].set_xlabel(r'$T_N$',fontsize=fs)

    #Save the figures
    plt.figure(fig1.number)
    plt.savefig(root_dir+alpha_dir+'ebtel2fl_L'+str(L)+'_tpulse'+str(t_pulse)+'_alpha'+str(alpha)+ '_' + species + '_heating_dem.eps',format='eps',dpi=1000)
    plt.figure(fig2.number)
    plt.savefig(root_dir+alpha_dir+'ebtel2fl_L'+str(L)+'_tpulse'+str(t_pulse)+'_alpha'+str(alpha)+ '_' + species + '_heating_TmaxVTn.eps',format='eps',dpi=1000)
    plt.figure(fig3.number)
    plt.savefig(root_dir+alpha_dir+'ebtel2fl_L'+str(L)+'_tpulse'+str(t_pulse)+'_alpha'+str(alpha)+ '_' + species + '_heating_hs_compare.eps',format='eps',dpi=1000)

    #Close the figures to avoid runtime warning of too many figures open
    plt.close('all')
