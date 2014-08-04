/***********************************************************************************

FILENAME: ebtel_functions_loop.c

AUTHOR Will Barnes

DATE: created: 7 March 2014

DESCRIPTION: This file contains all of the major functions used in the EBTEL model.
Additional functions can be found in ebtel_functions_util.c and ebtel_functions_param.c.
Full function descriptions can be found below. See ebtel_main.c for additional details.
***********************************************************************************/

//Include appropriate header file
#include "ebtel-2fl_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_loop_solver

FUNCTION DESCRIPTION: Given some heating input, this function solves the 0-D hydrodynamic equations spatially averaged over a loop half-length, computing the temperature, density, and pressure as functions of time. Given appropriate initial conditions (see below), the differential emission measure of the transition region and corona can also be calculated. 

INPUTS:
	ntot--total number of steps to be taken in the time integration
	loop_length--half-length of the coronal loop (cm) (top of chrmosphere to apex)
	total_time--total time over which the simulation will run (s)
   	opt--keyword structure: (set to ON (1) or OFF (0))
        dynamic--set to use dynamical r1 and r2
        dem_old--set to use old technique of computing DEM(T) in the TR
        rtv--set to use Rosner, Tucker, & Vaiana radiative loss function
        usage--this keyword is not optional and must be set to determine
               proper usage. See USAGE section in ebtel_main.c.
		solver--Euler solver (0) or RK solver (1)
		mode--IC calculation: static eq. (0), force with inputs (1), scaling laws (2)
		heating_shape--specifies which heating function to use
		index_dem--gives range over which the DEM is calculated
        Following value is not a keyword and contains actual values:
        energy_nt--mean energy of the nonthermal electrons (keV)
		T0--initial temperature (only used if mode=1) (K)
		n0--initial density (only used if mode=1) (cm^-3)
		h_nano--maximum value of the heating function (erg cm^-3 s^-1)
		t_pulse_half--time between onset of heating function and max heating (s)
		tau--time step (s)
OUTPUTS:
   	param_setter--pointer structure of arrays for all loop parameters; see ebtel_functions.h for a complete list of all structure members
***********************************************************************************/

struct ebtel_params_st *ebtel_loop_solver( int ntot, double loop_length, double total_time, struct Option *opt)
{
	/***********************************************************************************
								Variable Declarations
	***********************************************************************************/
	
	/*****Variable declarations*****/
	//Int
	int nk;
	int index_dem = opt->index_dem;
	int i;	//index over ntot
	int j;	//index over 451
	int k;	//index used for averaging over temporal DEM dimensions
	int flag_dem_tr;
	int j_min;
	int j_max;
	
	//Pointers
	double *kptr;
	double *state_ptr;
	double *log_tdem_ptr;
	double *ic_ptr;
	double *flux_ptr;
	
	//Double
	double r1;
	double r1_tr;
	double r12;
	double r12_tr;
	double r2;
	double r3;
	double r4;
	double r1e,r1i,r2e,r2i;
	
	double rad;
	double sc;
	
	double f_e;
	double f_i;
	double f_eq;
	double cf;
	
	double tau;
	double time = 0.;	//initialize time to zero
	double rad_loss;
	double t_e;
	double t_i;
	double n;
	double ta_e;
	double ta_i;
	double na;
	double p_e;
	double p_i;
	double pa_e;
	double pa_i;
	double v;
	
	double t_max;
	double t_min;
	double em;
	double dem0;
	double delta_t;
	double cond_e;
	double cond_i;
	double rad_cor;
	double rad_ratio;
	double f_ratio;
	
	//Array 
	double f_array[3];
	double state[15];
	double log_tdem[index_dem];
	double tdem[index_dem];
	double root_tdem[index_dem];
	double fourth_tdem[index_dem];
	double rad_dem[index_dem];
	double root_rad_dem[index_dem];
	double dem_cor_minus[ntot];
	double dem_tr_minus[ntot];
	double dem_tot_minus[ntot];
	
	//Two-dimensional array (dynamically allocate memory to avoid segmentation fault)
	double **dem_tr;
	dem_tr = malloc(ntot*sizeof(double *));
	double **dem_cor;
	dem_cor = malloc(ntot*sizeof(double *));
	for(i = 0; i<ntot; i++)
	{
		dem_tr[i] = malloc(index_dem*sizeof(double));
		dem_cor[i] = malloc(index_dem*sizeof(double));
	}
	
	//struct
	struct rk_params par;
	struct ebtel_rka_st *adapt;
	struct ebtel_params_st *param_setter = malloc(sizeof(struct ebtel_params_st));
	assert(param_setter != NULL);
	
	//Reserve memory for structure members that will be set
	param_setter->time = malloc(sizeof(double[ntot]));
	param_setter->tau = malloc(sizeof(double[ntot]));
	param_setter->heat = malloc(sizeof(double[ntot]));
	param_setter->temp_e = malloc(sizeof(double[ntot]));
	param_setter->temp_i = malloc(sizeof(double[ntot]));
	param_setter->ndens = malloc(sizeof(double[ntot]));
	param_setter->press_e = malloc(sizeof(double[ntot]));
	param_setter->press_i = malloc(sizeof(double[ntot]));
	param_setter->vel = malloc(sizeof(double[ntot]));
	param_setter->tapex_e = malloc(sizeof(double[ntot]));
	param_setter->tapex_i = malloc(sizeof(double[ntot]));
	param_setter->napex = malloc(sizeof(double[ntot]));
	param_setter->papex_e = malloc(sizeof(double[ntot]));
	param_setter->papex_i = malloc(sizeof(double[ntot]));
	param_setter->coeff_1 = malloc(sizeof(double[ntot]));
	param_setter->cond_e = malloc(sizeof(double[ntot]));
	param_setter->cond_i = malloc(sizeof(double[ntot]));
	param_setter->rad_cor = malloc(sizeof(double[ntot]));
	param_setter->rad = malloc(sizeof(double[ntot]));
		
	if(opt->usage == 4 || opt->usage == 1)
	{
		if(opt->usage == 4)
		{	
			param_setter->f_ratio = malloc(sizeof(double[ntot]));
			param_setter->rad_ratio = malloc(sizeof(double[ntot]));
		}
		param_setter->logtdem = malloc(sizeof(double[index_dem]));
		param_setter->dem_tr_log10mean = malloc(sizeof(double[index_dem]));
		param_setter->dem_cor_log10mean = malloc(sizeof(double[index_dem]));
		param_setter->dem_tot_log10mean = malloc(sizeof(double[index_dem]));
	}
	
	//DEBUG--reserve memory in param structure for parts of step in pe,pi
	param_setter->dpe = malloc(sizeof(double[ntot]));
	param_setter->dpe1 = malloc(sizeof(double[ntot]));
	param_setter->dpe2 = malloc(sizeof(double[ntot]));
	param_setter->dpe3 = malloc(sizeof(double[ntot]));
	param_setter->dpe4 = malloc(sizeof(double[ntot]));
	param_setter->dpe5 = malloc(sizeof(double[ntot]));
	param_setter->dpi = malloc(sizeof(double[ntot]));
	param_setter->dpi1 = malloc(sizeof(double[ntot]));
	param_setter->dpi2 = malloc(sizeof(double[ntot]));
	param_setter->dpi3 = malloc(sizeof(double[ntot]));
	
	/***********************************************************************************
									Initial Parameters
	***********************************************************************************/	
	
	//Calculate initial r values. In Klimchuk et al. (2008), these are denoted by c1,c2,c3. 
	//See ebtel_main.c documentation for correspondence between variables.
	r1 = ebtel_calc_c3();	//ratio of base to apex temperature; c3 in Klimchuk et al. (2008), Paper I from here on.
	r1_tr = r1;				//ratio of temperature at top of TR to apex temperature, can be set to 0.5 as well.
	r2 = ebtel_calc_c2();	//ratio of average to apex temperature; c2 in Paper I
	r3 = 2.0;				//ratio of TR to coronal radiative losses; c1 in Paper I
	r4 = 1.0;				//ratio of average to base velocity
	
	
	//Calculate ion and electron specific coefficients
	r1e = ebtel_calc_c3();
	r1i = ebtel_calc_c3();
	r2e = ebtel_calc_c2();
	r2i = ebtel_calc_c2();
	
	//Set temperature bins.
	//Lengths of the kpar array are different depending on the loss function we use.
	if (opt->rtv== 0)
	{nk = 7;
	}
	else
	{nk = 6;
	}
	//Declare kpar array which will set the bins based on our loss function choice
	double kpar[nk];
	//Call the ebtel_kpar_set function and pass the pointer kptr
	kptr = ebtel_kpar_set(opt->rtv);
	//Set the kpar array;
	for(i = 0; i < nk; i++)
	{
		kpar[i] = *(kptr + i);
	}
	//Calculate radiative loss function.
	rad = ebtel_rad_loss(1e+6,kpar,opt->rtv);
	
	/***********************************************************************************
						Set up DEM in Transition Region
	***********************************************************************************/
	if (opt->usage == 1 || opt->usage == 4)
	{	
		//Define temperature arrays for plotting and calculating DEM
		log_tdem_ptr = ebtel_linspace(0,index_dem-1,index_dem);		//set the linspace pointer using the ebtel_linspace function
		for(j = 0; j<index_dem; j++)
		{
			log_tdem[j] = *(log_tdem_ptr + j)/100 + 4;	//log of T in the region in which DEM is defined
			param_setter->logtdem[j] = log_tdem[j];		//Save logtdem to our parameter structure
			tdem[j] = pow(10,log_tdem[j]);				//T in the region in which DEM is defined
			root_tdem[j] = sqrt(tdem[j]);				//T^1/2 in the region in which DEM is defined
			fourth_tdem[j] =  pow(tdem[j],0.25);			//T^1/4 in the region in which DEM is defined
		}
		
		//These coefficients will be used in the old method of calculating DEM (global variables)
		ROOT_C2 = pow((KAPPA_0_E/(20*K_B)),0.5)/K_B;		//Calculate the root of c2 to avoid overflow in our calculation of dem_ev
		C3 = -5*K_B;
		C4 = pow((KAPPA_0_E/14),0.5)/K_B;
		
		//Radiation in the transition region. This loop just calculates the radiative loss function in the TR
		for(j = 0; j<index_dem; j++)
		{
			rad = ebtel_rad_loss(tdem[j],kpar,opt->rtv);		//Calculate the radiative loss function for temperature tdem[i]
			rad_dem[j] = rad;								//Set radiative loss function in the TR
			if (tdem[j] < 1e+4)
			{
				rad_dem[j] = 1;								//Check to see if we are outside the allowed temperature range
			}
			root_rad_dem[j] = sqrt(rad_dem[j]);
		}
	}
	
	/***********************************************************************************
							Initial Static Equilibrium
	***********************************************************************************/
	/*
	2 possible methods: (A) use EBTEL equilibrium (recommended) (B) use scaling laws (set by mode in opt structure)
	*/
	
	//Calculate initial temperature, density, pressure, and velocity.
	ic_ptr = ebtel_calc_ic(kpar,r3,loop_length,opt);
	r3 = *(ic_ptr + 0);
	rad = *(ic_ptr + 1);
	t_e = *(ic_ptr + 2);
	t_i = *(ic_ptr + 2);
	n = *(ic_ptr + 3);
	p_e = *(ic_ptr + 4);
	p_i = *(ic_ptr + 4);
	v = *(ic_ptr + 5);
	
	//Free the pointer used in the static equilibrium calculation
	free(ic_ptr);
	ic_ptr = NULL;
	
	//Set remaining initial parameters before iterating in time
	ta_e = t_e/r2e;
	ta_i = t_i/r2i;
	sc = ebtel_calc_lambda(t_e); //NOTE:Using both temperatures may not be right; na should be the same for both e,i since we assume ne = ni = n
	na = n*r2*exp(-2.0*loop_length/(PI*sc)*(1.0-sin(PI/5.0)));
	pa_e = K_B*na*ta_e;
	pa_i = KB_FACT*K_B*na*ta_i;
	
	//Set the initial values of our parameter structure
	param_setter->temp_e[0] = t_e;
	param_setter->temp_i[0] = t_i;
	param_setter->ndens[0] = n;
	param_setter->press_e[0] = p_e;
	param_setter->press_i[0] = p_i;
	param_setter->vel[0] = v;
	param_setter->tapex_e[0] = ta_e;
	param_setter->tapex_i[0] = ta_i;
	param_setter->napex[0] = na;
	param_setter->papex_e[0] = pa_e;
	param_setter->papex_i[0] = pa_i;
	param_setter->rad[0] = rad;
	param_setter->coeff_1[0] = r3;
	param_setter->heat[0] = ebtel_heating(time,opt);
	param_setter->time[0] = time;
	param_setter->tau[0] = opt->tau;
	
	//Print out the coefficients that we are starting the model with
	printf("********************************************************************\n");
	printf("Model Parameters\n");
	printf("Further modifications may be needed for ion and electron specific r1,r2 coefficients.\n");
	printf("r1 = %e\n",r1);
	printf("r1_e = %f\n",r1e);
	printf("r1_i = %f\n",r1i);
	printf("r2 = %e\n",r2);
	printf("r2_e = %f\n",r2e);
	printf("r2_i = %f\n",r2i);
	printf("r3 = %e\n",r3);
	printf("********************************************************************\n");
	
	//Print out the parameters that we are starting the model with
	printf("L = %e\n",loop_length);
	printf("Q = %e\n",param_setter->heat[0]);
	printf("T_e, T_i, Ta_e, Ta_i = %e, %e\n",param_setter->temp_e[0],param_setter->tapex_e[0]); //Apex e and ion temperatures may not always be equal at t=0
	printf("n, na = %e, %e\n",param_setter->ndens[0],param_setter->napex[0]);					//Depends on coefficients, equal for now
	printf("********************************************************************\n");
	
	/***********************************************************************************
							Time-dependent Heating
	***********************************************************************************/
	
	//Set the structure members of the par structure. This is used in both the RK and Euler methods.
	par.L = loop_length;
	par.kpar = kptr;
	par.r12 = r1/r2;
	par.r2 = r2;
	par.r4 = r4;
		
	//Set the initial timestep from opt structure.
	//This will be static if we are not using our adaptive solver
	tau = opt->tau;
	
	//Initialize the counter
	i = 0;
	
	//Begin the loop over the timesteps
	while(time < total_time && i < ntot)
	{	
		//Update the parameter structure
		par.q1 = ebtel_heating(time,opt);
		par.q2 = ebtel_heating(time+tau,opt);
		
		//Save the apex electron pressure to the par structure
		par.Pae = pa_e;
		
		//Update time
		time = time + tau;
		
		//Save time and heat to main data structure
		param_setter->heat[i+1] = par.q2;
		param_setter->time[i+1] = time;
		
		//DISABLE usage = 3 (NT e- flux ) for now
		//Set up non-thermal electron flux for usage option 3
		//if(opt.usage==3)
		//{
		//	par.flux_nt = loop_length*par.q2;
		//}
		//else
		//{
		//	par.flux_nt = 0;
		//}
		
		//Calculate radiative loss
		rad = ebtel_rad_loss(t_e,kpar,opt->rtv);
		param_setter->rad[i+1] = rad;

		//Calculate coefficients r1, r2, r3 (c3, c2, c1)
		r3 = ebtel_calc_c1(t_e,n,loop_length,rad);
		par.r3 = r3;
		param_setter->coeff_1[i+1] = r3;
		r2 = ebtel_calc_c2();
		par.r2 = r2;
		r1 = ebtel_calc_c3();
		r12 = r1/r2;
		par.r12 = r12;
		r12_tr = r1_tr/r2;
		
		//Unpack electron and ion heat flux
		flux_ptr = ebtel_calc_conduction(t_e,t_i,n,loop_length,rad,r3,opt->dynamic);
		f_e = *(flux_ptr + 0);
		f_i = *(flux_ptr + 1);
		f_eq = *(flux_ptr + 2);
		par.f_e = f_e;
		par.f_i = f_i;
		par.f_eq = f_eq;
		
 		//Free the flux pointer; it will be malloc'd on the next iteration
		free(flux_ptr);
		flux_ptr = NULL;
		
		/*****Step parameters forward in time using (1) RK method or (0) Euler stepper method*****/
		
		//Update the state vector
		state[0] = p_e;	
		state[1] = p_i;
		state[2] = n;
		state[3] = t_e;
		state[4] = t_i;
		
		if(opt->solver==0)	//Euler solver
		{	
			//Call the Euler routine
			state_ptr = ebtel_euler(state,tau,par);
		}
		else if(opt->solver==1)	//RK routine
		{	
			//Call the RK routine
			state_ptr = ebtel_rk(state,5,param_setter->time[i],tau,par,opt);	
		}
		else if(opt->solver==2)
		{	
			//Call the adaptive RK routine
			adapt = ebtel_rk_adapt(state,5,param_setter->time[i],tau,opt->error,par,opt);
			
			//Set the state vectore and timestep
			state_ptr = adapt->state;
			tau = adapt->tau;
		}
		
		//Calculate and store velocity
		v = ebtel_calc_vel(t_e,t_i,p_e,par);
		param_setter->vel[i+1] = v;		

		//Update p,n,t,tau and save to structure
		p_e = *(state_ptr + 0);
		param_setter->press_e[i+1] = p_e;
		p_i = *(state_ptr + 1);
		param_setter->press_i[i+1] = p_i;
		n = *(state_ptr + 2);
		param_setter->ndens[i+1] = n;
		t_e = *(state_ptr + 3);
		param_setter->temp_e[i+1] = t_e;
		t_i = *(state_ptr + 4);
		param_setter->temp_i[i+1] = t_i;
		
		//Save time step
		param_setter->tau[i+1] = tau;
		
		//DEBUG--save dpe,dpi steps
		param_setter->dpe[i] = *(state_ptr + 6);
		param_setter->dpe1[i] = *(state_ptr + 7);
		param_setter->dpe2[i] = *(state_ptr + 8);
		param_setter->dpe3[i] = *(state_ptr + 9);
		param_setter->dpe4[i] = *(state_ptr + 10);
		param_setter->dpe5[i] = *(state_ptr + 11);
		
		param_setter->dpi[i] = *(state_ptr + 12);
		param_setter->dpi1[i] = *(state_ptr + 13);
		param_setter->dpi2[i] = *(state_ptr + 14);
		param_setter->dpi3[i] = *(state_ptr + 15);
		
		//Free memory used by the state pointer. Free the adapt structure if we are using the adapt method.
		if(opt->solver==2)
		{
			//Free the structure as it will be malloc'd on the next go around
			free(adapt->state);
			adapt->state = NULL;
			free(adapt);
		}
		else
		{
			//Free the state ptr as it will be malloc'd again on the next go around
			free(state_ptr);
			state_ptr = NULL;
		}
		
		//Calculate new scale height
		sc = ebtel_calc_lambda(t_e); //NOTE: not completely sure about using sum of temperatures here
		
		//Calculate apex quantities
		ta_e = t_e/r2e;
		param_setter->tapex_e[i+1] = ta_e;
		ta_i = t_i/r2i;
		param_setter->tapex_i[i+1] = ta_i;
		na = n*r2*exp(-2.0*loop_length*(1.0-sin(PI/5.0))/(PI*sc));
		param_setter->napex[i+1] = na;
		pa_e = K_B*na*ta_e;
		param_setter->papex_e[i+1] = pa_e;
		pa_i = KB_FACT*K_B*na*ta_i;
		param_setter->papex_i[i+1] = pa_i;

		
		/*****Differential Emission Measure Calculation*****/
		//Check usage variable to determine whether we are calculating TR DEM
		if(opt->usage == 1 || opt->usage == 4)
		{	
			//Transition region
			if(r12_tr*t_e > tdem[index_dem-1])
			{
				printf(" Transition region T = %e K outside of DEM range\n",r12_tr*t_e);
				exit(0);
			}
			
			if(f_e != f_eq)
			{
				cf = f_e*f_eq/(f_e - f_eq);
			}
			else
			{
				cf = 1e+10*f_e;
			}
			
			//Make f array for ebtel_calc_tr_dem function
			f_array[0] = f_e;
			f_array[1] = f_eq;
			f_array[2] = cf;
			
			//Initialize dem_tr flag to zero. We want to know if dem_tr takes on negative values and if so we need to reset them.
			flag_dem_tr = 0;
			
			
			for(j=0; j<index_dem; j++)
			{
				//Check to see whether we are in the TR. If so, calculate dem_TR. Note: r12_tr*t[i] = T_0
				if( tdem[j] < r12_tr*t_e )
				{
						
					//Make call to function that calculates DEM for TR using method specified by opt.dem_old.
					dem_tr[i][j] = ebtel_calc_tr_dem(tdem[j],n,v,p_e,loop_length,sc,rad_dem[j],f_array,opt->dem_old);
					
					//Check whether the dem is less than zero and set the flag if the condition holds
					if(dem_tr[i][j] < 0)
					{
						flag_dem_tr = 1;
						printf("***** Negative DEM at i = %d\n", i);
						break;
					}
					
					//Set f_ratio value
					f_ratio = f_e/f_eq;
				}
			}
			
			//If we tripped the flag_dem_tr, then we need to reset dem_tr
			if(flag_dem_tr == 1)
			{
				for(j=0; j<index_dem; j++)
				{
					if(i != 0)
					{
						dem_tr[i][j] = dem_tr[i-1][j];
					}
					else
					{
						dem_tr[i][j] = 0;
					}
				}
			}
			
			//Corona (EM distributed uniformly over temperature interval [tmin,tmax])
			t_max = ebtel_max_val(t_e/r2e,1.1e+4);
			t_min = ebtel_max_val(t_e*(2.0 - 1/r2e),1e+4);
			
			j_max = (log10(t_max) - 4.0)*100;
			j_min = (log10(t_min) - 4.0)*100;
			
			//Calculate total emission measure
			em = 2*pow(n,2)*loop_length;			//factor of 2 for both legs
			
			delta_t = 1e4*(pow(10.0,((j_max + 0.5)/100.0)) - pow(10.0,((j_min - 0.5)/100.0)));
        	dem0 = em/delta_t;
        	
        	//Set every value between j_min and j_max DEM in the corona to dem0
        	for(j=j_min; j<=j_max; j++)
        	{
        		dem_cor[i][j] = dem0;
        	}
        	
        	//Transition region radiation losses based on DEM
        	if(opt->usage == 4)
        	{
        		rad_loss = 0;		//Initialize rad_loss to 0
        		
        		//Sum up radiative losses in the transition region
        		for(j=0; j<index_dem; j++)
        		{
        			if(tdem[j] < r12_tr*t_e)
        			{
        				rad_loss += dem_tr[i][j]*rad_dem[j]*tdem[j]*0.01*log(10.0);
        			}
        		}
        		
        		//Compute rad_ratio
        		rad_ratio = -rad_loss/f_eq;
        		param_setter->rad_ratio[i] = rad_ratio;
        		param_setter->f_ratio[i] = f_ratio;
        	}
		}
		
		//Set the conductive loss from the corona 
		cond_e = f_e;
		param_setter->cond_e[i] = cond_e;
		cond_i = f_i;
		param_setter->cond_i[i] = cond_i;
		//Set the coronal radiative loss value
		rad_cor = f_eq/r3;
		param_setter->rad_cor[i] = rad_cor;
		
		//Increment the counter
		i++;
	}
	
	//End of loop
	
	/***********************************************************************************
						Format Data and Free Pointers
	***********************************************************************************/
	
	//Set structure field to be used in file writer
	param_setter->i_max = i;
	
	//Format DEM data to be printed to file.
	//Take weighted time average for each T_DEM
	if(opt->usage == 1 || opt->usage == 4)
	{
		for(j = 0; j < index_dem; j++)
		{
			//Set the last entry that was left empty by the loop on t
			dem_tr[param_setter->i_max-1][j] = dem_tr[param_setter->i_max-2][j];
			dem_cor[param_setter->i_max-1][j] = dem_cor[param_setter->i_max-2][j];
			
			//Create a single dimensional array from a doubly indexed array
			for(k = 0; k<param_setter->i_max; k++)
			{
				dem_cor_minus[k] = dem_cor[k][j];
				dem_tr_minus[k] = dem_tr[k][j];
				dem_tot_minus[k] = dem_tr[k][j] + dem_cor[k][j];
			}
			
			//Compute the mean of each of our newly created arrays over the index i, compute the log10, and store it as an entry in the array dem_log10mean_{region}.
			param_setter->dem_cor_log10mean[j] = log10(ebtel_weighted_avg_val(dem_cor_minus,param_setter->i_max,param_setter->tau));
			param_setter->dem_tr_log10mean[j] = log10(ebtel_weighted_avg_val(dem_tr_minus,param_setter->i_max,param_setter->tau));
			param_setter->dem_tot_log10mean[j] = log10(ebtel_weighted_avg_val(dem_tot_minus,param_setter->i_max,param_setter->tau));
		}
	}
	
	//Free up memory used by ebtel_kpar_set and ebtel_linspace functions
	free(kptr);
	kptr = NULL;
	if(opt->usage==1 || opt->usage==4)
	{
		free(log_tdem_ptr);
		log_tdem_ptr = NULL;
	}
	for(i=0;i<ntot;i++)
	{
		free(dem_tr[i]);
		dem_tr[i] = NULL;
		free(dem_cor[i]);
		dem_cor[i] = NULL;
	}
	free(dem_tr);
	dem_tr = NULL;
	free(dem_cor);
	dem_cor = NULL;
	
	//DEBUG--set final values for debug parameters
	param_setter->dpe[param_setter->i_max] = param_setter->dpe[param_setter->i_max-1];
	param_setter->dpe1[param_setter->i_max] = param_setter->dpe1[param_setter->i_max-1];
	param_setter->dpe2[param_setter->i_max] = param_setter->dpe2[param_setter->i_max-1];
	param_setter->dpe3[param_setter->i_max] = param_setter->dpe3[param_setter->i_max-1];
	param_setter->dpe4[param_setter->i_max] = param_setter->dpe4[param_setter->i_max-1];
	param_setter->dpe5[param_setter->i_max] = param_setter->dpe5[param_setter->i_max-1];
	
	param_setter->dpi[param_setter->i_max] = param_setter->dpi[param_setter->i_max-1];
	param_setter->dpi1[param_setter->i_max] = param_setter->dpi1[param_setter->i_max-1];
	param_setter->dpi2[param_setter->i_max] = param_setter->dpi2[param_setter->i_max-1];
	param_setter->dpi3[param_setter->i_max] = param_setter->dpi3[param_setter->i_max-1];
	
	//Exit and return the structure that has been set appropriately
	return param_setter;
}

/***********************************************************************************

FUNCTION NAME: ebtel_kpar_set

FUNCTION DESCRIPTION: This function sets the kpar array to be used later in calculations of the radiative loss function. These are essentially the temperature bins used in calculating the radiative loss function. Either the RK or RTV method is used.

INPUTS:
	rtv_opt--use either the RK or RTV method 
OUTPUTS:
	*kpar--pointer for the kpar array

***********************************************************************************/

double * ebtel_kpar_set(int rtv_opt)
{	
	//Check option input to decide which method to use
	if (rtv_opt == 0)
	{	
		double *kpar = malloc(sizeof(double[7]));
		//Raymond-Klimchuk Loss function
		kpar[0] = 1.0e+4;
		kpar[1] = 9.3325e4;
		kpar[2] = 4.67735e5;
        kpar[3] = 1.51356e6;
        kpar[4] = 3.54813e6;
        kpar[5] = 7.94328e6;
        kpar[6] = 4.28048e7;
        
        return kpar;
	}
	else
	{
		double *kpar = malloc(sizeof(double[6]));
		//Rosner-Tucker-Vaiana Loss function
		kpar[0] = pow(10.0,4.3);
        kpar[1] = pow(10.0,4.6);
        kpar[2] = pow(10.0,4.9);
        kpar[3] = pow(10.0,5.4);
        kpar[4] = pow(10.0,5.75);
        kpar[5] = pow(10.0,6.3);
        
        return kpar;
	}
}

/***********************************************************************************

FUNCTION NAME: ebtel_rad_loss

FUNCTION_DESCRIPTION: This function calculates the radiative loss function using either the Raymond-Klimchuk (RK) method or the Rosner-Tucker-Vaiana (RTV) method. Adopted from pro radLoss from the original EBTEL IDL code.

INPUTS:
	temp--temperature (K)
	kpar--holds the kpar structure if it has been previously set; placeholder otherwise.
	rtv_opt--option to use the RTV method (1) or the RK method (1).
	
OUTPUTS:
	rad--radiative loss function

***********************************************************************************/

double ebtel_rad_loss( double temp, double kpar[], int rtv_opt)
{
	//Declare rad 
	double rad;
	
	//Find out which method is being used
	if (rtv_opt == 0)
	{
		//RK loss function
    	if ( temp > kpar[6] ){ 
        	rad = 1.96e-27*pow(temp,0.5);
        }
    	else if ( temp > kpar[5] ){ 
        	rad = 5.4883e-16/temp;
        }
    	else if ( temp > kpar[4] ){
        	rad = 3.4629e-25*(pow(temp,1./3.));
        }
    	else if ( temp > kpar[3] ){ 
        	rad = 3.5300e-13/(pow(temp,1.5));
        }
    	else if ( temp > kpar[2] ){ 
        	rad = 1.8957e-22;
        }
    	else if ( temp > kpar[1] ){
        	rad = 8.8669e-17/temp;
        }
    	else if ( temp >= kpar[0] ){
        	rad = 1.0909e-31*pow(temp,2.);
        }
    	else{
        	rad = 0.0;
        }
	}
	else
	{
		//RTV loss function
		if (temp > kpar[5]){ 
        	rad = pow(10.,-17.73)/pow(temp,2./3.);
        }
    	else if (temp > kpar[4] ){
        	rad = pow(10.,-21.94);
        }
    	else if (temp > kpar[3] ){
        	rad = pow(10.,-10.4)/pow(temp,2.);
        }
    	else if (temp > kpar[2] ){ 
        	rad = pow(10.,-21.2);
        }
    	else if ( temp > kpar[1] ){
        	rad = pow(10.,-31.0)*pow(temp,2.);
        }
    	else if ( temp >= kpar[0] ){ 
        	rad = pow(10.,-21.85);
        }
    	else{
        	rad = 0.0;
        }
	}
	
	return rad;
}



/***********************************************************************************

FUNCTION NAME: ebtel_calc_tr_dem

FUNCTION_DESCRIPTION: This function calculates the differential emission measure in the transition region using the new method described in Paper II. The method used here is to rewrite the energy equation for the transition region in terms of dtds such thatit is quadratic in dtds. The parameters aaa,bbb,ccc are the coefficients on dtds in this quadratic equation. 

INPUTS:
	tdem--temperature in the TR for which we calculate the DEM
	n--density at time t_i
	v--velocity at time t_i
	p--pressure at time t_i
	L--loop half length
	sc--scale height at temperature T_i
	rad_dem--radiative loss in the transition region
	f_a--array holding f,f_eq,cf
	option--option to either use old (1) or new method (0)
	
OUTPUTS:
	dem_tr--differential emission measure in the transition region

***********************************************************************************/

double ebtel_calc_tr_dem(double tdem, double n, double v, double p, double L, double sc, double rad_dem, double f_a[3], int option)
{
	double dem_tr;	//Declare what will be returned by the function

	//First check to see what method we are using to calculate the DEM in the TR
	if(option==0)
	{
		/*********New Method*********/
		//Declare variables
		double a;
		double b;
		double c;
		double p2kt2;
		double dtds1;
		double dtds2;
		double dtds;
	
		//Calculate necessary coefficients for quadratic formula
		a = KAPPA_0_E*pow(tdem,1.5);
		b = -5.0*K_B*n*v;
		p2kt2 = pow((p*exp(2.0*sin(PI/5.0)*L/(PI*sc))/(2.0*K_B*tdem)),2.0);	//(p/2kT)^2; calculate this here for convenience
		c = -p2kt2*rad_dem;
		dtds1 = (-b + sqrt(pow(b,2.0) - 4.0*a*c))/(2.0*a);
		dtds2 = (-b - sqrt(pow(b,2.0) - 4.0*a*c))/(2.0*a);
		dtds = ebtel_max_val(dtds1,dtds2);
		dem_tr = 2.0*p2kt2/dtds;		//factor of 2 for both legs of the loop
	}
	else
	{
		/*********Old Method*********/
		//Declare variables
		double dem_ev;
		double dem_con;
		double dem_eq;
		double root_tdem = sqrt(tdem);
		double root_rad_dem = sqrt(rad_dem);
		double fourth_tdem = pow(tdem,0.25);
		
		//Approximation to TR reg. DEM when evaporation dominates
		dem_ev = (ROOT_C2/n)*(ROOT_C2*pow(p,2)/root_tdem)/v;
		//Approximation to TR reg. DEM when condensation dominates
		dem_con = C3*n*v/rad_dem;
		//Approximation to TR reg. DEM under equilibrium conditions
		dem_eq = C4*p/(root_rad_dem*fourth_tdem);

		//Calculate DEM for the transition region
		dem_tr = 2.*(f_a[0]*dem_ev + f_a[2]*dem_eq - f_a[1]*dem_con)/(f_a[0] + f_a[2] - f_a[1]); //factor of 2 for both legs
	}
	
	return dem_tr;

}
 