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
 
 Function name: ebtel_euler
 
 Function description: This function implements a Euler in EBTEL. It uses a simple 
 Euler stepper method to solve our simplified hydrostatic equations.
  
 Input
	s--vector that stores the current state of the system
	tau--timestep
	par--structure that holds necessary parameters
	opt--structure that holds necessary input parameters
 
 Return
 	s_out--updated state vector
	
 *********************************************************************************/
 
 double * ebtel_euler(double s[], double tau, struct rk_params par, struct Option opt)
 {
 	//Declare variables
 	double p_e;
	double p_i;
 	double n,n_old;
 	double T_e;
	double T_i;
 	double dn;
 	double dp_e;
	double dp_i;
	double p_ev;
	double nu_ei;
	double vdPds_TR;
	double vdPds_C;
 	double *s_out = malloc(sizeof(double[5]));
 
 	//Unravel the state vector
	//p_e and n are set to old value so that we are consistent at which time t we are evaluating our expressions
 	p_e = s[0];
	p_i = s[1];
 	n_old = s[2];
	T_e = s[3];
 	T_i = s[4];
	
	//Calculate enthalpy flux
	p_ev = 2./3.*(par.f_eq - par.f_e);
	
	//Calculate collisional frequency
	nu_ei = ebtel_collision_freq(T_e,T_i,n_old);
	
	//Approximate TR and C integrals of v*dPe/ds terms
	vdPds_TR = 0;
	vdPds_C = par.v*par.Pae; 
 
	//Advance n in time
	//NOTE: At this point we have not changed the coefficients r1, r2, r3 so these expressions may change 
	dn = (p_ev/(par.L*K_B*par.r12*T_e))*tau;
	n = n_old + dn;
	
	//Advance p_e,p_i in time
	dp_e = (2./3.*(par.q1 + 1./par.L*par.f_eq*(1. + 1./par.r3) + 1./par.L*(vdPds_TR + vdPds_C)) + K_B*n_old*nu_ei*(T_i - T_e))*tau;
	p_e = p_e + dp_e;
	
	dp_i = (2./3./par.L*(vdPds_TR + vdPds_C) + K_B*n_old*nu_ei*(T_e - T_i))*tau;
	p_i = p_i + dp_i;
	
	//Calculate T
	T_e = p_e/(n*K_B);
	T_i = p_i/(n*K_B);
	
	//Update the state vector and return it
	s_out[0] = p_e;
	s_out[1] = p_i;
	s_out[2] = n;
	s_out[3] = T_e;
	s_out[4] = T_i;
	
	return s_out;

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
 
 double * ebtel_rk(double s[], int n, double t, double tau, struct rk_params par,struct Option opt)
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
 	f1 = ebtel_rk_derivs(s,t,0,par,opt);
	//Make the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f1 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}
 	
 	//Compute the second function f2
 	f2 = ebtel_rk_derivs(s_temp,t_half,1,par,opt);
	//Rebuild the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f2 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}
 	
 	//Compute the third function f3
 	f3 = ebtel_rk_derivs(s_temp,t_half,1,par,opt);
	//Rebuild the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f3 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}
 	
 	//Compute the fourth function f4
 	f4 = ebtel_rk_derivs(s_temp,t_full,2,par,opt);
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
 	free(f2);
 	free(f3);
 	free(f4);
 	
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
 
 struct ebtel_rka_st *ebtel_rk_adapt(double s[], int n, double t, double tau, double err, struct rk_params par, struct Option opt)
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
 			scale = err*(fabs(s_small_2[j]) + fabs(s_big[j]))/2.0;
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
 
/**********************************************************************************
 
 Function name: ebtel_rk_derivs
 
 Function description: This function solves the derivatives for the EBTEL model. More
 specifically, it computes dpdt, dndt, dTdt at some given time t and returns these 
 derivatives to the ebtel_rk function so that the RK routine can be executed.
  
 Input
	s--vector that stores the current state of the system
	t--current time
	tau_opt--tells what kind of time step is being used to appropriately compute heating
			 (0)t: q1
			 (1)t + half_tau: (q1 + q2)/2
			 (2)t + tau: q2
	par--structure that holds necessary parameters
 	opt--structre that holds input parameters
 
 Return
 	derivs--updated derivative state vector
	
 *********************************************************************************/
 
 double * ebtel_rk_derivs(double s[], double t, int tau_opt, struct rk_params par, struct Option opt)
 {
	
 	//Declare variables
 	double p_e,p_i;
 	double n;
	double v;
 	double T_e,T_i;
 	double rad;
 	double r3;
 	double f_e,f_i,f_eq;
	double nu_ei;
 	double q;
	double p_ev;
	double vdPds_TR,vdPds_C;
 	double dp_edt;
	double dp_idt;
 	double dndt;
 	double dT_edt;
	double dT_idt;
 	double *derivs = malloc(sizeof(double[5]));
 	int nk;
 	int i;
	
	double *flux_ptr;
 
 	//Unravel the state vector
	//p_e and n are set to old value so that we are consistent at which time t we are evaluating our expressions
 	p_e = s[0];
	p_i = s[1];
 	n = s[2];
	T_e = s[3];
 	T_i = s[4];
 	
 	//Make the kpar array
 	if(opt.rtv==0)
 	{
 		nk = 7;
 	}
 	else
 	{
 		nk = 6;
 	}
 	double kpar[nk];
 	for(i=0; i<nk; i++)
 	{
 		kpar[i] = *(par.kpar + i);
 	}
 	
 	//Compute the radiative loss function 
 	rad = ebtel_rad_loss(T_e,kpar,opt.rtv);
 	
 	//Compute the coefficient r3
 	r3 = ebtel_calc_c1(T_e,n,par.L,rad);
 	
 	//Compute heat flux
	flux_ptr = ebtel_calc_conduction(T_e,T_i,n,par.L,rad,r3,opt.dynamic);
	f_e = *(flux_ptr + 0);
	f_i = *(flux_ptr + 1);
	f_eq = *(flux_ptr + 2);
	free(flux_ptr);
	flux_ptr = NULL;
	
	//Set the heating depending on the time input. For tau_half, we average between the i and i+1 heating
	if(tau_opt==0)
	{
		q = par.q1;
	}
	else if(tau_opt==1)
	{
		q = (par.q1 + par.q2)/2;
	}
	else
	{
		q = par.q2;
	}
	
	//Calculate the enthalpy flux
	p_ev = 2./3.*(f_eq - f_e);
	
	//Calculate velocity
	v = p_ev/p_e*par.r4;
	
	//Approximate TR and C integrals of v*dPe/ds terms
	vdPds_TR = p_ev;
	vdPds_C = v*p_e; 
	
	//Calculate collision frequency
	nu_ei = ebtel_collision_freq(T_e,T_i,n);
	
	//Now compute the derivatives of each of the quantities in our state vector
	dp_edt = (2./3.*(q + 1./par.L*f_eq*(1. + 1./r3) + 1./par.L*(vdPds_TR + vdPds_C)) + K_B*n*nu_ei*(T_i - T_e));
	dp_idt = (2./3./par.L*(vdPds_TR + vdPds_C) + K_B*n*nu_ei*(T_e - T_i));
	
	dndt = (p_ev/(par.L*K_B*par.r12*T_e));
	
	dT_edt = T_e*(1./p_e*dp_edt - 1./n*dndt);
	dT_idt = T_i*(1./p_i*dp_idt - 1./n*dndt);
	
	//Set the derivative state vector
	derivs[0] = dp_edt;
	derivs[1] = dp_idt;
	derivs[2] = dndt;
	derivs[3] = dT_edt;
	derivs[4] = dT_idt;	
	
	//Return the pointer
	return derivs;
 	
 }