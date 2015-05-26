#ebtel2fl_configure_main.py

#Will Barnes
#11 May 2015

#Import needed modules
import argparse
import sys
sys.path.append('../../../bin/')
from ebtel2fl_configure import Configurer

#Declare parser object
parser = argparse.ArgumentParser(description='Script that prints configuration files for EBTEL-2fluid runs.')
#Add arguments to parser
parser.add_argument("-s","--species",help="Species to which the heating was applied for particular run.")
parser.add_argument("-as","--amp_switch",help="Switch to decide between power-law and uniform heating.")
parser.add_argument("-a","--alpha",type=float,help="Spectral index for the power-law distribution used.")
parser.add_argument("-L","--loop_length",type=float,help="Loop half-length.")
parser.add_argument("-t","--t_pulse",type=float,help="Width of the heating pulse used for the particular run.")
parser.add_argument("-S","--solver",help="Solver used to compute solutions.")
#Declare the parser dictionary
args = parser.parse_args()

#Set heating parameters
Q0 = 1e+23 #lower bound on nanoflare energy
Q1 = 1e+25 #upper bound on nanoflare energy
Ah = 1e+14 #loop cross sectional area
Hn = 8.3e-3 #Average nanoflare energy distributed over the total time

#Configure all static dictionary options
config_dict = {'usage_option':'dem','rad_option':'rk','dem_option':'new','heat_flux_option':'limited','solver':args.solver,'ic_mode':'st_eq'}
config_dict['total_time'] = 80000
config_dict['tau'] = 1.0
config_dict['rka_error'] = 1.0e-6
config_dict['index_dem'] = 451
config_dict['sat_limit'] = 0.166667
config_dict['h_back'] = 3.4e-6
config_dict['heating_shape'] = 'triangle'
config_dict['t_start_switch'] = 'file'
config_dict['t_end_switch'] = 'file'
config_dict['T0'] = 1.0e+6
config_dict['n0'] = 1.0e+8
config_dict['t_start'] = 0.0
config_dict['t_pulse_half'] = 0.5*args.t_pulse
config_dict['mean_t_start'] = 1000
config_dict['std_t_start'] = 1000

#Configure directory-level parameters
config_dict['heat_species'] = args.species
config_dict['amp_switch'] = args.amp_switch
config_dict['alpha'] = args.alpha
config_dict['loop_length'] = args.loop_length
config_dict['amp0'] = Q0/(config_dict['loop_length']*1.0e+8*Ah*args.t_pulse) #lower bound on nanoflare volumetric heating rate
config_dict['amp1'] = Q1/(config_dict['loop_length']*1.0e+8*Ah*args.t_pulse) #upper bound on nanoflare volumetric heating rate

config = Configurer(config_dict,'/data/datadrive2/EBTEL-2fluid_runs/',Hn=Hn,mc=5000)
config.vary_wait_time(250,5000,250)
config.print_job_array_config()
