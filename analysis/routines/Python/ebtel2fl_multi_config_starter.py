#Name: ebtel2fl_multi_config_starter.py
#Author: Will Barnes
#Date: 21 February 2015

#Description: Loop through sets of directories to start multiple runs of the EBTEL 2-fluid model

#Set the base directories
base_dir_config = '/data/datadrive2/EBTEL-2fluid_runs/ion_heating_runs/'
base_dir_exec = '/home/wtb2/Documents/EBTEL-2fluid_repo/'

#Import the necessary modules
import sys
sys.path.append(base_dir_exec + 'bin/')
import ebtel2fl_wrapper as ew

#Create the tuple that will vary our directory
var_tuple = (1.5,2.0,2.5,'uniform')

#Loop over the tuple that varies the directories and run over all config files in every directory
for i in range(len(var_tuple)):
    #Create the directory name
    var_dir = 'alpha' + str(var_tuple[i]) + '/config/'
    #Run all config files in this directory
    ew.run_ebtel(base_dir_exec + 'bin/',base_dir_config+var_dir)
    #Print the status to the screen
    print "Status: ",var_dir," runs completed."
