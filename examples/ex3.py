"""
Third example: single-fluid, single asymmetric pulse
"""

import sys
import os
import subprocess

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
sns.set_context('notebook',font_scale=1.5)

top_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(top_dir,'rsp_toolkit/python'))
from xml_io import InputHandler,OutputHandler

def run_example():

    #configure the run
    ih = InputHandler(os.path.join(top_dir,'config','ebtel.example.cfg.xml'))
    base_dir = ih.lookup_vars()
    base_dir['calculate_dem'] = True
    base_dir['output_filename'] = os.path.join(top_dir,'examples','ex3')
    base_dir['heating']['partition'] = 0.5
    base_dir['force_single_fluid'] = True
    base_dir['heating']['events'] = [{'event':{'rise_start':0.0,'rise_end':100.0,'decay_start':500.0,'decay_end':1500.0,'magnitude':0.005}}]
    #print the file
    oh = OutputHandler(base_dir['output_filename']+'.xml',base_dir)
    oh.print_to_xml()

    #run the model
    subprocess.call([os.path.join(top_dir,'bin','ebtel++.run'),'-c',base_dir['output_filename']+'.xml'])

    #bring in data and preprocess dem results
    results = np.loadtxt(base_dir['output_filename'])
    results_dem_tr = np.loadtxt(base_dir['output_filename']+'.dem_tr')
    results_dem_corona = np.loadtxt(base_dir['output_filename']+'.dem_corona')
    results_dem_temperature = results_dem_tr[0,:]
    results_dem_total = np.average(results_dem_tr[1:,:]+results_dem_corona[1:,:],axis=0,weights=np.gradient(results[:,0]))
    results_dem_tr = np.average(results_dem_tr[1:,:],axis=0,weights=np.gradient(results[:,0]))
    results_dem_corona = np.average(results_dem_corona[1:,:],axis=0,weights=np.gradient(results[:,0]))

    #setup figure
    fig = plt.figure(figsize=(9,6))
    gs = gridspec.GridSpec(3,2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0],sharex=ax1)
    ax3 = fig.add_subplot(gs[2,0],sharex=ax1)
    ax4 = fig.add_subplot(gs[:,1])

    #plasma parameters
    ax1.plot(results[:,0],results[:,-1])
    ax2.plot(results[:,0],results[:,1]/1e+6,color=sns.color_palette()[0],label=r'$T_e$')
    ax2.plot(results[:,0],results[:,2]/1e+6,color=sns.color_palette()[2],label=r'$T_i$')
    ax3.plot(results[:,0],results[:,3]/1e+8)
    #dem
    ax4.plot(results_dem_temperature,results_dem_tr,color=sns.color_palette()[0],label=r'$\mathrm{DEM}_{TR}$')
    ax4.plot(results_dem_temperature,results_dem_corona,color=sns.color_palette()[2],label=r'$\mathrm{DEM}_{C}$')
    ax4.plot(results_dem_temperature,results_dem_total,color=sns.color_palette()[1],label=r'$\mathrm{DEM}_{total}$')

    #axes options
    ax1.set_ylabel(r'$H$ (erg cm$^{-3}$ s$^{-1}$)')
    ax1.tick_params(axis='x',labelbottom='off')
    ax1.locator_params(axis='y',tight=True,nbins=5)
    ax2.set_ylabel(r'$T$ (MK)')
    ax2.tick_params(axis='x',labelbottom='off')
    ax2.locator_params(axis='y',tight=True,nbins=5)
    ax3.set_ylabel(r'$n$ ($10^8$ cm$^{-3}$)')
    ax3.set_xlabel(r'$t$ (s)')
    ax3.locator_params(axis='y',tight=True,nbins=5)
    ax3.set_xlim([results[0,0],results[-1,0]])
    ax2.legend(loc='best')
    ax4.set_xlabel(r'$T$ (K)')
    ax4.set_ylabel(r'$\mathrm{DEM}$ (cm$^{-5}$ K$^{-1}$)')
    ax4.set_xlim([10**(4.5),10**(7.5)])
    ax4.set_ylim([10**(20.0),10**(23.5)])
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.legend(loc='best')

    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)),'ex3.png'))
if __name__ == '__main__':
    run_example()
