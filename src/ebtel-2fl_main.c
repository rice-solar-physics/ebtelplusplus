/************************************************************************************

FILENAME: ebtel-2fl_main.c

AUTHOR Will Barnes

DATE: created: 30 June 2014

DESCRIPTION: This file makes the function call to the primary function in ebtel-2fl_functions_loop.c.
It first defines some basic parameters of the EBTEL model including the input parameters as well
as the constants used throughout. Following the necessary function call, it saves the data.

-------------------------------------------------------------------------------------
EBTEL (Enthalpy Based Thermal Evolution of Loops) computes 0-D
hydrodynamic equations.  This software was originally developed by Klimchuk et al,. (2008) and
later improved upon by Cargill et al., (2012a,b). It solves the set of hydrostatic equations by
averaging over the loop half-length and then integrating in time. Additional details are available
in the two references given above.

Note on variable correspondence with Klimchuk et al. (2008)
   r1 = c_3
   r2 = c_2
   r3 = c_1
   f, ff = F_0
   f_eq, ff_eq = - R_tr
   dem_eq = DEM_se

INTENSITIES:
   For observations in which temperature response function, G(T), has units of
   DN s^-1 pix^-1 cm^5 and the loop diameter, d, is larger than the pixel dimension, l_pix:

      I_cor_perp = d/(2L)* Int{G(T)*dem_cor(T)*dT}
      I_tr_perp = d/l_pix * Int{G(T)*dem_tr(T)*dT}
      I_tr_parallel = Int{G(T)*dem_tr(T)*dT} ,

   for lines-of-sight perpendicular and parallel to the loop axis.  I_tr_perp assumes that
   the transition region is thinner than l_pix.

MISCELLANEOUS COMMENTS:
   Runs much more quickly if the transition region DEM is not computed.
   Speed can be increased by increasing the minimum DEM temperature from 10^4 to, say, 10^5 K
      or by decreasing the maximum DEM temperature from 10^8.5 to, say, 10^7.5
      (search on 450 and 451).
   The equilibrium base heat flux coefficient of 2/7 is appropriate for uniform heating;
      a coefficient of 4/7 is more appropriate for apex heating.
   To have equal amounts of thermal and nonthermal heating:  flux_nt = heat*length.
   It is desirable to have a low-level background heating during the cooling phase so that the
      coronal temperature does not drop below values at which the corona DEM is invalid.
   r1 = c_3 = 0.7 gives more accurate coronal evolution than the original 0.5, especially in
      the late phase of cooling.  However, it produces excess DEM at the very hottest temperatures
      during impulsive events, and the transition region DEM is somewhat elevated.  We have
      therefore introduced r1_tr = 0.5, which provides a more accurate transition region DEM at
      the same time that r1 = 0.7 provides a more accurate radiative cooling.
   v = (c_3/c_2)*(t_tr/t)*v_0 = (r1/r2)*(t_tr/t)*(v/r4) at temperature t_tr in the transition
      region, where t is the average coronal temperature.

USAGE:
(set in usage field of opt parameter structure. See input list above.
(1)--include transition region DEM
     additional outputs: dem_tr, dem_cor, logtdem
(2)--exclude transition region DEM (faster)
(3)--include nonthermal electron energy flux
     additional inputs: flux_nt, energy_nt
(4)--compute rad_ratio (25% more computing time)
     additional outputs: dem_tr, dem_cor, logtdem, f_ratio, rad_ratio

************************************************************************************/

#include "ebtel-2fl_functions.h"

int main (int argc, char *argv[])
{
	//Use clock to time the entire EBTEL program
	clock_t time_start;
	clock_t time_diff;
	double time_elapsed;

	//Start the timer
	time_start = clock();

	/************************************************************************************
								Variable Declarations
	************************************************************************************/

	/******Variable Declarations******/
	//Struct
	struct ebtel_params_st *params_final;		//Declare instance of structure ebtel_params_st
	struct Option *opt;

	//Global definitions (declarations in ebtel_functions.h)
	KAPPA_0_E = 7.8e-7;			//Spitzer coefficient for electron thermal conduction
	KAPPA_0_I = 3.2e-8;			//Spitzer coefficient for ion thermal conduction
	M_EL = 9.11e-28;			//mass of e- in grams
	Q_E = 4.8032e-10;			//charge of e- in stat coloumbs
	G0_SUN = 2.74e+4;			//gravitational acceleration at the solar surface in cm/s^2
	PI = 3.14159265359;
	TWO_SEVENTHS = 2./7.;
	SEVEN_HALVES = 3.5;
	TWO_THIRDS = 2./3.;
	GAMMA = 5.0/3.0;

	//Set global variables based on He/H abundance
	ebtel_calc_abundance();

	int i,n;
	int quiet_flag = 0;
	double L;
	char filename_in[250];
	double *kptr;

	/**********************************
	Read configuration file
	**********************************/

	//Read in parameters from file
	//Set default filename
	sprintf(filename_in,"../config/ebtel_config.xml");
	//Check if a filename was specified at the command line
	for(i = 0; i<argc; i++)
	{
		//Read in filename
		if(strcmp(argv[i],"quiet")!=0 && i>0)
		{
			sprintf(filename_in,"%s",argv[i]);
		}
		else if(strcmp(argv[i],"quiet")==0)
		{
			//Raise the quiet flag
			quiet_flag = 1;
		}
	}

	//Pass the filename to the strcut setter to read inputs into opt structure
	opt = ebtel_input_setter(filename_in);

	/************************************************************************************
									Initial Parameters
	************************************************************************************/

	//Set total number of steps using the initial timestep and total time
	if(strcmp(opt->solver,"euler")==0 || strcmp(opt->solver,"rk4")==0)
	{
		//For static timesteps, this is just the total time divided by the timestep
		n = ceil(opt->total_time/opt->tau) + 1;
	}
	else if(strcmp(opt->solver,"rka4")==0)
	{
		//When using the adaptive method, this is a guess and additional memory will be allocated if necessary
 	   	n = ceil(opt->total_time) + 1;
	}
	else
	{
		printf("Invalid solver option.\n");
		exit(0);
	}

	//Define loop half-length and change to appropriate units
	L = 1.0e+8*opt->loop_length;	//convert from Mm to cm

	//Set members of the Option opt structure
	opt->energy_nt = 8.01e-8;	//50 keV in ergs
	
	//Set temperature bins for calculating radiative loss function
 	//Make the kpar array
	NK=6;	//There is a correspondence with the length of KPAR[]; if NK ever changes, change length of KPAR[] in header file
	kptr = ebtel_kpar_set(opt->rad_option);
 	for(i=0; i<NK; i++)
 	{
 		KPAR[i] = *(kptr + i);
 	}
	free(kptr);
	kptr=NULL;

	/************************************************************************************
									Heating
	************************************************************************************/

	//Configure start times, end times and amplitudes for selected heating events. These
	//arrays can be configured using the parameters specified in the above input file, using
	//normal and power-law distributions, or through additional input files.
	ebtel_heating_config(opt,filename_in);

	/************************************************************************************
									Start the Model
	************************************************************************************/

	//Print a header to the screen that gives the user input information
	if(quiet_flag == 0)
	{
		ebtel_print_header(n, opt);
	}

	//Make the call to the ebtel_loop_solver function. This function sets the members of the structure params_final. Each member
	//is a pointer to an array
	params_final = ebtel_loop_solver(n, L, opt);

	/************************************************************************************
									Save the Data
	************************************************************************************/

	//Write the contents of params_final to a file. See output for filename.
	ebtel_file_writer(opt, params_final);

	//Count the number of events and print it to the screen
	int num_q_events = ebtel_count_events(params_final,opt);
	printf("Number of simulated heating events: %d\n",num_q_events);

	/****************Done writing data to file. Free up memory reserved by pointers.******************/

	//Free up memory used by the structure params_final and Option
	ebtel_free_mem(params_final,opt);

	//Stop the timer
	time_diff = clock() - time_start;
	time_elapsed = time_diff*1000/CLOCKS_PER_SEC;

	//Time elapsed
	printf("The process took %f milliseconds to run\n",time_elapsed);

	//Exit with no errors
	return 0;

}
