#ebtel2fl_run_main.py

#Will Barnes
#11 May 2015

#Import needed modules
import argparse
import numpy as np
import sys
sys.path.append('/home/wtb2/Documents/EBTEL-2fluid_repo/bin/')
from ebtel2fl_run import Runner

#Declare parser object
parser = argparse.ArgumentParser(description='Script to perform Monte Carlo-like EBTEL simulations by running over same parameter set multiple times')

#Add arguments to parser
parser.add_argument("-k","--job_index",type=int,help="Index for (Tn,i) pair to specify waiting time and run number from job array.")
parser.add_argument("-s","--species",help="Species to which the heating was applied for particular run.")
parser.add_argument("-a","--alpha",help="Spectral index for the power-law distribution used. Use absolute value in this instance.")
parser.add_argument("-L","--loop_length",help="Loop half-length.")
parser.add_argument("-t","--t_pulse",help="Width of the heating pulse used for the particular run.")
parser.add_argument("-S","--solver",help="Solver used to compute solutions.")

#Declare the parser dictionary
args = parser.parse_args()

#set executable directory
exec_dir = '/home/wtb2/Documents/EBTEL-2fluid_repo/bin/'
#set top level directory
top_dir = '/data/datadrive2/EBTEL-2fluid_runs/'
top_dir = top_dir + args.species+'_heating_runs/alpha'+str(args.alpha)+'/config/'

#Load job array file
job_array = np.loadtxt(top_dir+'ebtel2fl_L'+str(args.loop_length)+'_tpulse'+str(args.t_pulse)+'_'+args.solver+'_job_array.conf',dtype='int')[args.job_index]

#Instantiate class
sub_dir = 'ebtel2fl_L'+str(args.loop_length)+'_tn%d_tpulse'+str(args.t_pulse)+'_'+args.solver
eb_run = Runner(exec_dir,top_dir+sub_dir%(job_array[0])+'/')
eb_run.run_ebtel_single(sub_dir%(job_array[0])+'_'+str(job_array[1])+'.xml')