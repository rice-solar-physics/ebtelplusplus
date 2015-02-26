#Name: ebtel2fl_configure_runs
#Author: Will Barnes
#Date: 17 February 2014

#Description: Configure EBTEL-2fluid runs

#Set up root directory extension
root = '/data/datadrive2/EBTEL-2fluid_runs/'
root_ebtel2fl = '/home/wtb2/Documents/EBTEL-2fluid_repo/'

#Import necessary modules
import sys
sys.path.append(root_ebtel2fl + 'bin/')
import ebtel2fl_wrapper as ew
import numpy as np

#Define the function that configures the start, end time arrays
def config_start_end_time(t_wait, t_total, t_pulse):
    #Calculate the number of pulses
    N = int(np.ceil(t_total/(t_pulse + t_wait)))

    #Create the needed arrays
    t_start_array = np.empty([N])
    t_end_array = np.empty([N])

    #Create start and end time arrays
    for i in range(N):
        t_start_array[i] = i*(t_pulse + t_wait)
        t_end_array[i] = t_start_array[i] + t_pulse

    return {'num_events':N,'t_start_array':t_start_array,'t_end_array':t_end_array}

def power_law_dist(x0, x1, x, alpha):
    return ((x1**(alpha+1) - x0**(alpha+1))*x + x0**(alpha+1))**(1/(alpha+1))


#Set heating parameters
Q0 = 1e+23 #lower bound on nanoflare energy
Q1 = 1e+25 #upper bound on nanoflare energy
Ah = 1e+14 #loop cross sectional area
Hn = 8.3e-3 #Average nanoflare energy distributed over the total time

#Set up array of wait times
T_wait = np.arange(250,5250,250)
#Set the pulse duration
t_pulse = 200.0

#Configure static parameters
config_dict = {'usage_option':'dem','rad_option':'rk','dem_option':'new','heat_flux_option':'limited','solver':'euler','ic_mode':'st_eq'}
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
config_dict['t_pulse_half'] = 50.0
config_dict['mean_t_start'] = 1000
config_dict['std_t_start'] = 1000

#Configure directory-level parameters
config_dict['heat_species'] = 'ion'
config_dict['amp_switch'] = 'uniform'
config_dict['alpha'] = -1.5
config_dict['loop_length'] = 20.0
config_dict['amp0'] = Q0/(config_dict['loop_length']*1.0e+8*Ah*t_pulse) #lower bound on nanoflare volumetric heating rate
config_dict['amp1'] = Q1/(config_dict['loop_length']*1.0e+8*Ah*t_pulse) #upper bound on nanoflare volumetric heating rate

#Set up directory to print config files
top_dir = config_dict['heat_species']+'_heating_runs/'
if config_dict['amp_switch'] == 'uniform':
    top_dir = top_dir + 'alpha' + config_dict['amp_switch'] + '/'
else:
    top_dir = top_dir + 'alpha' + str(-1*config_dict['alpha']) + '/'
config_dir = root + top_dir + 'config/'
data_dir = root + top_dir + 'data/'

#Loop over different values of time between successive nanoflares
for i in range(len(T_wait)):

    #Calculate the start and end time arrays and the number of events
    heat_times = config_start_end_time(T_wait[i], config_dict['total_time'], t_pulse)
    config_dict['num_events'] = heat_times['num_events']
    config_dict['start_time_array'] = heat_times['t_start_array']
    config_dict['end_time_array'] = heat_times['t_end_array']

    #Set the uniform peak nanoflare energy (for triangular pulses)
    config_dict['h_nano'] = 2*Hn*config_dict['total_time']/(config_dict['num_events']*t_pulse)

    #Concatenate the filename
    fn = 'ebtel2fl_L' + str(config_dict['loop_length']) + '_tn' + str(T_wait[i]) + '_tpulse' + str(t_pulse) + '_' + config_dict['solver']

    #Set the ouput filename
    config_dict['output_file'] = data_dir + fn

    #Print the config file (use same filename as output with _config)
    ew.print_xml_config(config_dict,config_file=config_dir+fn+'.xml')
