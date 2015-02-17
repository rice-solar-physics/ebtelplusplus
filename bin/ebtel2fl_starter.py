#/opt/local/bin/python
#ebtel_starter.py

#Uses the ebtel_wrapper module to configure an EBTEL-C run and plot the results

import ebtel2fl_wrapper as ew
import numpy as np

#Heating file extension
root = '/Users/willbarnes/Documents/Rice/Research/EBTEL-2fluid_repo/'
heat_ext = '/Users/willbarnes/Documents/Rice/Research/EBTEL_repo/analysis/data/'

#Create dictionary with desired parameters
#Switches
run_dictionary = {'usage_option':'dem','rad_option':'rk','dem_option':'new','heat_flux_option':'limited','solver':'rka4','ic_mode':'st_eq','heat_species':'electron'}
#General input
run_dictionary['total_time'] = 10000
run_dictionary['tau'] = 1.0
run_dictionary['loop_length'] = 75.0
run_dictionary['rka_error'] = 1.0e-6
run_dictionary['index_dem'] = 451
run_dictionary['T0'] = 0.85e+6
run_dictionary['n0'] = 0.36e+8
run_dictionary['sat_limit'] = 1.0/6.0
#Heating parameters
run_dictionary['heating_shape'] = 'triangle'
run_dictionary['t_start_switch'] = 'uniform'
run_dictionary['t_end_switch'] = 'uniform'
run_dictionary['amp_switch'] = 'uniform'
run_dictionary['num_events'] = 1
run_dictionary['t_start'] = 0.0
run_dictionary['t_pulse_half'] = 250.0
run_dictionary['h_nano'] = 1.5e-3
run_dictionary['h_back'] = 3.4e-6
run_dictionary['mean_t_start'] = 15000
run_dictionary['std_t_start'] = 20000
run_dictionary['alpha'] = -2
run_dictionary['amp0'] = 0.001
run_dictionary['amp1'] = 0.01
if run_dictionary['t_start_switch']=='file':
    run_dictionary['start_time_array'] = np.loadtxt(heat_ext+'hydrad_warren_start.txt')
if run_dictionary['t_end_switch']=='file':
    run_dictionary['end_time_array'] = np.loadtxt(heat_ext+'hydrad_warren_end.txt')
if run_dictionary['amp_switch']=='file':
    run_dictionary['amp_array'] = np.loadtxt(heat_ext+'hydrad_warren_amp.txt')

#Specify config filename
config_file = 'ebtel-2fl_config_run1.xml'

#Construct data filename based on inputs
data_file =  'ebtel-2fldatL' + str(run_dictionary['loop_length']) + '_' + run_dictionary['usage_option'] + '_' + run_dictionary['heating_shape'] + '_' + run_dictionary['solver'] + '.txt'
data_file_dem =  'ebtel-2fldemdatL' + str(run_dictionary['loop_length']) + '_' + run_dictionary['usage_option'] + '_' + run_dictionary['heating_shape'] + '_' + run_dictionary['solver'] + '.txt'

#Construct figure filename
fig_name = '/Users/willbarnes/Documents/Rice/Research/EBTEL-2fluid_figures/ion_heating_runs/ebtel-2fl-C_'+run_dictionary['heat_flux_option']+'_'+run_dictionary['solver']+'_tau'+str(run_dictionary['tau'])+'_'+run_dictionary['heat_species'] + '_' + run_dictionary['heating_shape'] + '_N' + str(run_dictionary['num_events']) + '_' + run_dictionary['t_start_switch'] + '.eps'

#Print the configuration file
ew.print_xml_config(run_dictionary,config_file=root+'config/'+config_file)

#Run EBTEL-C executable
ew.run_ebtel(root+'bin/','../config/',config_file=config_file)

#Plot the results
ew.plot_ebtel(root+'data/',data_file)#,print_fig_filename=fig_name)
ew.plot_ebtel_dem(root+'data/',data_file_dem)

