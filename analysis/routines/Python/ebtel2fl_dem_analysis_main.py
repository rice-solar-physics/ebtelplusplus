#ebtel2fl_dem_analysis_main.py

#Will Barnes
#15 May 2015

#Import needed modules
import sys
import argparse
sys.path.append('../../../bin/')
import ebtel2fl_dem as ebd
import ebtel2fl_plot_em as ebpe

#set root directory
root_dir = '/data/datadrive2/EBTEL-2fluid_runs/'
#set figure file format
figname = root_dir+'%s_heating_runs/alpha%f/ebtel2fl_L%f_tpulse%f_alpha%f_%s_heating'

#make parameter vectors
loop_length = [20.0,40.0,60.0,80.0,100.0,120.0]
alpha = [1.5,2.0,2.5]

#set static parameters
tpulse = 100.0
solver = 'rka4'
mc = 100

#parse species argument
parser = argparse.ArgumentParser(description='Script that performs DEM analysis for EBTEL-2fluid runs.')
#Add arguments to parser
parser.add_argument("-s","--species",help="Species to which the heating was applied for particular run.")
#Declare the parser dictionary
args = parser.parse_args()

#iterate over variable parameters
for i in range(len(alpha)):
    for j in  range(len(loop_length)):
        dema = ebd.DEMAnalyzer(root_dir,args.species,alpha[i],loop_length[j],tpulse,solver,mc=mc)
        dema.process_raw()
        dema.many_slopes()
        demp = ebpe.DEMPlotter(dema.temp_em,dema.em,alpha[i])
        demp.plot_em_curves(print_fig_filename=figname%(args.species,alpha[i],loop_length[j],tpulse,alpha[i],args.species)+'_dem')
        demp.plot_em_max(print_fig_filename=figname%(args.species,alpha[i],loop_length[j],tpulse,alpha[i],args.species)+'_TmaxVTn')
        demp.plot_em_slopes(dema.a_cool,dema.a_hot,print_fig_filename=figname%(args.species,alpha[i],loop_length[j],tpulse,alpha[i],args.species)+'_hs_compare')