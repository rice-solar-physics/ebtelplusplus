/***********************************************************************************

FILENAME: ebtel_functions_solvers.c

AUTHOR Will Barnes

DATE: created: 19 April 2014

DESCRIPTION: This file contains functions used for the Euler and fourth-order Runge-Kutta
solvers implemented in EBTEL. An additional routine is included for the adaptive step
option that can be chosen in ebtel_main.

***********************************************************************************/

//Include appropriate header file
#include "ebtel-2fl_functions.h"

/**********************************************************************************

 Function name: ebtel_derivs

 Function description: This function solves the derivatives for the EBTEL model. More
 specifically, it computes dpdt, dndt, dTdt at some given time t and returns these
 derivatives to the ebtel_rk function so that the RK routine can be executed. Additionally, if the Euler
 solver option is specified, it will return the updated state vector p_e,p_i,n,T_e,T_i.

 Input
	s--vector that stores the current state of the system
	t--current time
	par--structure that holds necessary parameters
 	opt--structre that holds input parameters

 Return
 	derivs--updated derivative state vector (or state vector for Euler solver)

 *********************************************************************************/

 double * ebtel_derivs(double s[], double t, struct rk_params par, struct Option *opt)
 {
	//Declare variables
	double p_e,p_i;
	double n;
	double T_e,T_i;
	double rad;
	double r1e,r1i,r2e,r2i;
	double r3;
	double xi;
	double f_e,f_i,f_eq;
	double nu_ei;
	double qi,qe;
	double p_ev;
	double R_tr;
	double vdPds_TR;
	double dp_edt;
	double dp_idt;
	double dndt;
	double dT_edt;
	double dT_idt;
	double *derivs = malloc(sizeof(double[6]));
	double *flux_ptr;

	//Unpack state vector
	p_e = s[0];
	p_i = s[1];
	n = s[2];
	T_e = s[3];
	T_i = s[4];

	//Compute the radiative loss function
	rad = ebtel_rad_loss(T_e,opt->rad_option);

	//Compute the coefficient r3
	r3 = ebtel_calc_c1(T_e,T_i,n,par.L,rad,opt);

	//Compute heat flux
	flux_ptr = ebtel_calc_conduction(T_e,T_i,n,par.L,rad,r3,opt->sat_limit,opt->heat_flux_option);
	f_e = *(flux_ptr + 0);
	f_i = *(flux_ptr + 1);
	f_eq = *(flux_ptr + 2);
	free(flux_ptr);
	flux_ptr = NULL;

	//Set the heating--check whether this is ion or electron heating
	if(strcmp(opt->heat_species,"electron")==0)
	{
		qe = ebtel_heating(t,opt);
		qi = 0.0;
	}
	else if(strcmp(opt->heat_species,"ion")==0)
	{
		qi = ebtel_heating(t,opt);
		qe = 0.0;
	}
	else
	{
		printf("Invalid heat species option.\n");
		exit(0);
	}

	//Calculate collisional frequency
	nu_ei = ebtel_collision_freq(T_e,T_i,n);

	//Calculate ratio of base temperatures for ions and electrons
	r2e = ebtel_calc_c2();
	r2i = ebtel_calc_c2();
	r1e = ebtel_calc_c3();
	r1i = ebtel_calc_c3();
	xi = r1e/r1i*r2i/r2e*T_e/T_i/KB_FACT;

	//Calculate the radiative loss of the transition region
	R_tr = -f_eq;

	//Approximate TR integral of v*dPe/ds term
	vdPds_TR = (f_e - xi*f_i + R_tr)/(1. + xi);

	//Calculate enthalpy flux
	p_ev = (GAMMA - 1.)/GAMMA*(vdPds_TR - f_e - R_tr);

	//Advance n in time
	dndt = (r2e/(K_B*par.L*r1e*T_e)*p_ev);

	//Advance p_e,p_i in time
	dp_edt = (GAMMA - 1.)*(qe - 1./par.L*R_tr*(1. + 1./r3) + 1./par.L*vdPds_TR) + K_B*n*nu_ei*(T_i - T_e);
	dp_idt = (GAMMA - 1.)*(qi - 1./par.L*vdPds_TR) + KB_FACT*K_B*n*nu_ei*(T_e - T_i);

	dT_edt = T_e*(1./p_e*dp_edt - 1./n*dndt);
	dT_idt = T_i*(1./p_i*dp_idt - 1./n*dndt);

	//Return updated parameters (Euler) or derivatives (RK)
	if(strcmp(opt->solver,"euler")==0)
	{
		derivs[0] = dp_edt*(opt->tau) + p_e;
		derivs[1] = dp_idt*(opt->tau) + p_i;
		derivs[2] = dndt*(opt->tau)+ n;
		derivs[3] = derivs[0]/(derivs[2]*K_B);
		derivs[4] = derivs[1]/(derivs[2]*KB_FACT*K_B);
	}
	else if(strcmp(opt->solver,"rk4")==0 || strcmp(opt->solver,"rka4")==0)
	{
		derivs[0] = dp_edt;
		derivs[1] = dp_idt;
		derivs[2] = dndt;
		derivs[3] = dT_edt;
		derivs[4] = dT_idt;
	}
	else
	{
		printf("Invalid solver option.\n");
		exit(0);
	}

	return derivs;
 }

  /**********************************************************************************

 Function name: ebtel_rk

 Function description: This function implements a Runge-Kutta routine in EBTEL. It will
 call the ebtel_rk_derivs function to compute derivatives of the necessary functions.
 This implementation is based on a routine written in MATLAB by A. L. Garcia (see
 Numerical Methods for Physics, Prentice Hall, 1994).

 Input
	s--vector that stores the current state of the system
	n--number of variables in vector s
	t--current time
	tau--timestep
	par--structure that holds necessary parameters
	opt--structure that contains necessary input parameters

 Return
 	s_out--updated state vector

 *********************************************************************************/

 double * ebtel_rk(double s[], int n, double t, double tau, struct rk_params par,struct Option *opt)
 {
 	//Declare variables
 	double half_tau;
 	double t_half;
 	double t_full;
 	double f_temp;
 	double *f1;
 	double *f2;
 	double *f3;
 	double *f4;
 	double s_temp[n];
 	double *s_out = malloc(sizeof(double[n]));
 	int i;

 	//Set time variables
 	half_tau = 0.5*tau;
 	t_half = t + half_tau;
 	t_full = t + tau;

 	//Compute the first function f1
 	f1 = ebtel_derivs(s,t,par,opt);
	//Make the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f1 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}

 	//Compute the second function f2
 	f2 = ebtel_derivs(s_temp,t_half,par,opt);
	//Rebuild the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f2 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}

 	//Compute the third function f3
 	f3 = ebtel_derivs(s_temp,t_half,par,opt);
	//Rebuild the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f3 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}

 	//Compute the fourth function f4
 	f4 = ebtel_derivs(s_temp,t_full,par,opt);
	//Rebuild the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f4 + i);
 		s_temp[i] = s[i] + tau*f_temp;
 	}

 	//Now compute the final resulting state vector
 	for(i=0; i<n; i++)
 	{
 		s_out[i] = s[i] + 1./6.*tau*(*(f1+i) + *(f4+i) + 2.*( *(f2+i) + *(f3+i) ));
 	}

 	//Free memory of all f functions
 	free(f1);
	f1 = NULL;
 	free(f2);
	f2 = NULL;
 	free(f3);
	f3 = NULL;
 	free(f4);
	f4 = NULL;

 	//Return the final resulting state vector (pointer)
 	return s_out;

 }

 /**********************************************************************************

 Function name: ebtel_rk_adapt

 Function description: This function implements an adaptive step size for the Runge-Kutta
 routine used in EBTEL. It returns a structure pointer containing the current time, timestep,
 and state vector. It calls the ebtel_rk function. This implementation is based on a routine
 written in MATLAB by A. L. Garcia (see Numerical Methods for Physics, Prentice Hall, 1994).

 Input
 	s--state vector
 	n--size of s
 	t--current time
 	tau--current time step
 	err--desired truncation error
 	par--structure used in computing derivatives
 	opt--structure containing input options

 Return
	rka_params--structure containing updated time, state, and timestep

*********************************************************************************/

 struct ebtel_rka_st *ebtel_rk_adapt(double s[], int n, double t, double tau, struct rk_params par, struct Option *opt)
 {
 	/**Declare variables**/
 	//Int
 	int i;
 	int j;
 	int max_try = 100;

 	//Double
 	double safe1 = 0.9;
 	double safe2 = 1.1;
	double safe3 = 4;
 	double scale;
 	double x_diff;
 	double error_ratio;
 	double t_save = t;
 	double time;
 	double half_tau;
 	double old_tau;
	double tau_tc;
 	double epsilon = 1.0e-16;

	//Pointers
 	double *x_small_1;
 	double *x_small_2;
 	double *x_big;

	//Arrays
 	double s_small_1[n];
 	double s_small_2[n];
 	double s_big[n];

	//Structure
 	struct ebtel_rka_st *rka_params = malloc(sizeof(struct ebtel_rka_st));
 	//Reserve memory for structure members that will be set
 	rka_params->state=malloc(sizeof(double[n]));

 	//Loop over maximum number of attempts to satisfy error bound
 	for(i=0; i<max_try; i++)
 	{
 		//Take two small steps

 		//First small step
 		half_tau = 0.5*tau;

		x_small_1 = ebtel_rk(s,n,t_save,half_tau,par,opt);

 		//Unpack the x_small_1 pointer
 		for(j=0;j<n;j++)
 		{
 			s_small_1[j] = *(x_small_1 + j);
 		}

 		//Update the time vector
 		time = t_save + half_tau;

 		//Second small step
		x_small_2 = ebtel_rk(s_small_1,n,time,half_tau,par,opt);

		//Unpack the x_small_2 pointer
 		for(j=0;j<n;j++)
 		{
 			s_small_2[j] = *(x_small_2 + j);
 		}

 		//Take single big step
 		x_big = ebtel_rk(s,n,t_save,tau,par,opt);

		//Unpack the x_big pointer
 		for(j=0;j<n;j++)
 		{
 			s_big[j] = *(x_big + j);
 		}

 		//Update the time vector
 		time = t_save + tau;

		error_ratio = 0.;
		//Compute estimated truncation error
		for(j=0;j<n;j++)
 		{
			scale = opt->rka_error*(fabs(s_small_2[j]) + fabs(s_big[j]))/2.0;
			x_diff = s_small_2[j] - s_big[j];
 			//Return the maximum value of the error ratio
			error_ratio = ebtel_max_val(error_ratio,fabs(x_diff)/(scale + epsilon));
		}

 		//Estimate new tau value (including safety factors)
 		old_tau = tau;
 		tau = safe1*old_tau*pow(error_ratio,-1./5.);
 		tau = ebtel_max_val(tau,old_tau/safe2);

 		//If error is acceptable, set our structure values and return the structure
 		if(error_ratio < 1)// && error_ratio != 0)
 		{
 			//Set the structure fields
			tau = ebtel_min_val(tau,safe3*old_tau);

			//Check thermal conduction timescale
			tau_tc = ebtel_thermal_conduction_timescale(s[3],s[4],s[2],opt);

			tau = ebtel_min_val(tau,0.5*tau_tc);

 			rka_params->tau = tau;
 			for(j=0;j<n;j++)
 			{
 				rka_params->state[j] = s_small_2[j];
 			}

 			//Free the memory of the pointers that were used
 			free(x_small_1);
 			x_small_1 = NULL;
 			free(x_small_2);
 			x_small_2 = NULL;
 			free(x_big);
 			x_big = NULL;

 			//Return the structure
 			return rka_params;
 		}

 		//Free memory of the pointers that were used. They will be malloc'd on the next iteration
 		free(x_small_1);
 		x_small_1 = NULL;
 		free(x_small_2);
 		x_small_2 = NULL;
 		free(x_big);
 		x_big = NULL;
 	}

 	//If we've finished the loop without meeting the error requirement, then return an error
 	printf("Error: Adaptive Runge-Kutta routine failed after %d iterations. Exiting the program\n",i);

	//Exit if the routine fails
	exit(0);

 }
