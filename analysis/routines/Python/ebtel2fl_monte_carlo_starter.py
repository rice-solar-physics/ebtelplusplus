#Name: ebtel2fl_monte_carlo_starter.py
#Author: Will Barnes
#Date: 21 February 2015
#Modified: 27 March 2015

#Description: Execute N_mc EBTEL-2fluid runs for a single point in the (L,alpha,species,Tn) parameter space

#Set the base directories
base_dir = '/data/datadrive2/EBTEL-2fluid_runs/'
base_dir_exec = '/home/wtb2/Documents/EBTEL-2fluid_repo/'

#Import the necessary modules
import sys
sys.path.append(base_dir_exec + 'bin/')
import ebtel2fl_wrapper as ew
import argparse

#Set up the parser to get output from command line
#Declare parser object
parser = argparse.ArgumentParser(description='Script to perform Monte Carlo-like EBTEL simulations by running over same parameter set multiple times')

#Add arguments to parser
parser.add_argument("-s","--species",help="Species to which the heating was applied for particular run.")
parser.add_argument("-a","--alpha",type=float,help="Spectral index for the power-law distribution used. Use absolute value in this instance")
parser.add_argument("-L","--loop_length",type=float,help="Loop half-length.")
parser.add_argument("-t","--t_pulse",type=float,help="Width of the heating pulse used for the particular run.")
parser.add_argument("-S","--solver",help="Solver used to compute solutions.")
parser.add_argument("-T","--t_wait",help="Waiting time between consecutive nanoflares")

#Declare the parser dictionary
args = parser.parse_args()

#Set the Monte-Carlo number--number of simulations per parameter space coordinate
N_mc = 100

#Preallocate space for the DEM arrays
dem_temp = ew.np.empty(451,5)

#Set the config directory
var_dir = args.species + '_heating_runs/alpha' + str(args.alpha) + '/'

#Set the config file
file_prefix = 'ebtel2fl_L'+str(args.loop_length)+'_tn'+str(args.t_wait)+'_tpulse'+str(args.t_pulse)+'_'+args.solver

#Loop over the tuple that varies the directories and run over all config files in every directory
for i in range(N_mc):
    #Run the executable for the given config file
    ew.run_ebtel(base_dir_exec + 'bin/',base_dir+var_dir+'config/',config_file=file_prefix+'.xml')
    #Load the DEM file for this run and add it to the total
    dem_temp = dem_temp + ew.np.loadtxt(base_dir+var_dir+'data/'+file_prefix+'_dem.txt')

#Average the N_mc measurements and print them to a file
dem_temp = dem_temp/float(N_mc)
ew.np.savetxt(base_dir+var_dir+'data/'+file_prefix+'_dem.txt',dem_temp)

#Print the status to the screen
print "MC Status: species = " + args.species +", L = " + str(args.loop_length) + ", alpha = "+str(args.alpha)+",Tn = "+str(args.t_wait) +" runs completed."
