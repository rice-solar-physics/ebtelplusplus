#/opt/local/bin/python
#ebtel_starter.py

#Uses the ebtel_wrapper module to configure an EBTEL-C run and plot the results

import ebtel_wrapper as ew
import numpy as np

#Heating file extension
root = '/Users/willbarnes/Documents/Rice/Research/EBTEL_repo/'
heat_ext = root+'analysis/data/'

#Create dictionary with desired parameters
#Switches
run_dictionary = {'usage_option':'tr','rad_option':'rtv','dem_option':'new','heat_flux_option':'classical','solver':'rka4','ic_mode':'force'}
#General input
run_dictionary['total_time'] = 2000
run_dictionary['tau'] = 2.0
run_dictionary['loop_length'] = 25
run_dictionary['rka_error'] = 1.0e-6
run_dictionary['index_dem'] = 451
run_dictionary['T0'] = 7.36e+5
run_dictionary['n0'] = 1.27e+8
#Heating parameters
run_dictionary['heating_shape'] = 'square'
run_dictionary['t_start_switch'] = 'uniform'
run_dictionary['t_end_switch'] = 'uniform'
run_dictionary['amp_switch'] = 'power_law'
run_dictionary['num_events'] = 10
run_dictionary['t_start'] = 10.0
run_dictionary['t_pulse_half'] = 50.0
run_dictionary['h_nano'] = 0.01
run_dictionary['h_back'] = 1.1e-6
run_dictionary['mean_t_start'] = 5000
run_dictionary['std_t_start'] = 2000
run_dictionary['alpha'] = -2
run_dictionary['amp0'] = 0.001
run_dictionary['amp1'] = 0.01
run_dictionary['start_time_array'] = np.loadtxt(heat_ext+'hydrad_warren_start.txt')
run_dictionary['end_time_array'] = np.loadtxt(heat_ext+'hydrad_warren_end.txt')
run_dictionary['amp_array'] = np.loadtxt(heat_ext+'hydrad_warren_amp.txt')

#Specify config filename
config_file = 'ebtel_config_test.xml'
#Construct data filename based on inputs
data_file =  'ebteldatL' + str(run_dictionary['loop_length']) + '_' + run_dictionary['usage_option'] + '_' + run_dictionary['heating_shape'] + '_' + run_dictionary['solver'] + '.txt'

#Print the configuration file
ew.print_xml_config(run_dictionary,config_file=root+'config/'+config_file)

#Run EBTEL-C executable
ew.run_ebtel(root+'bin/','../config/',config_file=config_file)

#Plot the results
ew.plot_ebtel(root+'data/',data_file)

