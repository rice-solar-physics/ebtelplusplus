#ebtel2fl_run_main.py

#Will Barnes
#11 May 2015

#Import needed modules
import argparse
from ebtel2fl_run import Runner

#Declare parser object
parser = argparse.ArgumentParser(description='Script to perform Monte Carlo-like EBTEL simulations by running over same parameter set multiple times')

#Add arguments to parser
parser.add_argument("-s","--species",help="Species to which the heating was applied for particular run.")
parser.add_argument("-a","--alpha",help="Spectral index for the power-law distribution used. Use absolute value in this instance")
parser.add_argument("-L","--loop_length",type=float,help="Loop half-length.")
parser.add_argument("-t","--t_pulse",type=float,help="Width of the heating pulse used for the particular run.")
parser.add_argument("-S","--solver",help="Solver used to compute solutions.")

#Declare the parser dictionary
args = parser.parse_args()

#set executable directory
exec_dir = '/Users/willbarnes/Documents/Rice/Research/EBTEL-2fluid_repo/bin/'
#set top level directory
top_dir = '/Users/willbarnes/Documents/Rice/Research/EBTEL-2fluid_repo/analysis/routines/Python/'
top_dir = top_dir + args.species+'_heating_runs/alpha'+str(args.alpha)+'/config/'

#Instantiate class
eb_run = Runner(exec_dir,top_dir)
sub_dir = 'ebtel2fl_L'+str(args.loop_length)+'_tn'+str(250)+'_tpulse'+str(args.t_pulse)+'_'+args.solver+'/'
eb_run.run_ebtel_multi_parallel(sub_dir=sub_dir)