/***********************************************************************************

FILENAME: ebtel_functions_param.c

AUTHOR Will Barnes

DATE: created: 7 March 2014

DESCRIPTION: This file contains all of the functions used to calculate important parameters
(c1,c2,c3,lambda) in EBTEL. All functions are described below. For additional details,
see ebtel_main.c.

***********************************************************************************/

//Include appropriate header file
#include "ebtel-2fl_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_calc_c1

FUNCTION_DESCRIPTION: This function calculates the ratio between radiative losses from the transition region and the corona The correspondence between the EBTEL code and Klimchuk et al.(2008) is c1=r3. Adapted from the function pro cacl_c1 in the original IDL EBTEL code.

INPUTS:
	t--current time (s)
	temp_e--electron temperature (K)
	temp_i--ion temperature (K)
	den--electron number density (cm^-3).
	llength--loop half length (cm).
	rad--radiative loss function
	opt--structure that holds all input options

OUTPUTS:
	c1--ratio of transition region to coronal radiative loss functions

***********************************************************************************/

double  ebtel_calc_c1(double temp_e, double temp_i, double den, double loop_length, double rad, struct Option *opt )
{

	//Declare variables
	double sc;
	double n_eq_2,noneq2;
	double r2,r3;
	double grav_correction = 1.0;
	double loss_correction = 1.0;
	double r3_eqm_0 = 2.0;		//value in equilibrium with no gravity, -2/3 loss power law; hardcoded, should NOT be user-specified
	double l_fact = 5.0;		//geometric factors for inclusion of gravitational effects, l_fact^-1 = s/L*1/2 where <s/L> approx 0.4 so 1/l_fact apprrox 0.2 or 1/5
	
	//Calculate the scale height
	sc = ebtel_calc_lambda(temp_e,temp_i);
	
	//Calculate r2 value
	r2 = ebtel_calc_c2();

	
	//Adjust values for gravity
	if (strcmp(opt->r3_grav_correction,"true")==0 || strcmp(opt->r3_grav_correction,"True")==0)
	{
		grav_correction = exp(4*sin(PI/l_fact)*loop_length/(PI*sc));
	}
	
	//Adjust for loss function
	if (strcmp(opt->r3_loss_correction,"true")==0 || strcmp(opt->r3_loss_correction,"True")==0)
	{
		loss_correction = 1.95e-18/pow(temp_e,TWO_THIRDS)/rad;
	}
	
	//Calculate over/under density
	n_eq_2 = KAPPA_0_E*pow((temp_e/r2),SEVEN_HALVES)/(SEVEN_HALVES*r3_eqm_0*grav_correction*loss_correction*rad*pow(loop_length,2));
	noneq2 = pow(den,2)/n_eq_2;
	
	//Use different values of r3 based on value of noneq2
	if (noneq2 < 1.0)
	{
		//Hot loops equilibrium value of c1
		r3 = (2.0*r3_eqm_0 + opt->r3_cond_0*(1.0/noneq2 - 1.0))/(1.0 + 1.0/noneq2);
	}
	else
	{	
		//Radiative loops transition from equilibrium (noneq2-->1)
		r3 = (2.0*r3_eqm_0 + opt->r3_rad_0*(noneq2 - 1.0))/(1.0 + noneq2);
	}
	
	//Return the value of the parameter
	return r3*grav_correction*loss_correction;
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_c2

FUNCTION_DESCRIPTION: This function calculates the parameter c2 which is the ratio of average temperature to apex temperature. Adapted from the function pro calc_c2 in the original IDL EBTEL code.

INPUTS:

OUTPUTS:
	c2--ratio of tbar to ta

***********************************************************************************/

double ebtel_calc_c2( void )
{
	double c2;

	c2 = 0.9;
	//c2 = 0.87;	//Paper I value

	return c2;
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_c3

FUNCTION_DESCRIPTION: This function calculates the parameter c3 which is the ratio of base temperature to apex temperature. Adapted from the function pro calc_c3 in the original IDL EBTEL code.

INPUTS:

OUTPUTS:
	c3--ratio of t0 to ta

***********************************************************************************/

double ebtel_calc_c3(void)
{
	double c3;

	c3 = 0.6;
	//c3 = 0.5;	//Paper I value

	return c3;

}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_sound_speed

FUNCTION_DESCRIPTION: This function calculates the sound speed in the coronal plasma.

INPUTS:
	temp_e--electron temperature (K)
	temp_i--ion temperature (K)
	
OUTPUTS:
	cs--sound speed (cm s^-1)

***********************************************************************************/

double ebtel_calc_sound_speed(double temp_e, double temp_i)
{
	return pow(5.0/3.0*KB_FACT*K_B*(temp_e+temp_i)/(M_P),0.5);
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_lambda

FUNCTION_DESCRIPTION: This function calculates the scale height for a given temperature. Adapted from the function pro calc_lambda in the original IDL EBTEL code.

INPUTS:
	temp_e--electron temperature (K)
	temp_i--ion temperature (K)

OUTPUTS:
	sc--scale height (cm)

***********************************************************************************/

double ebtel_calc_lambda( double temp_e, double temp_i )
{
	double sc;

	sc = (KB_FACT*K_B*(temp_e + temp_i)/M_P)/G0_SUN;

	return sc;
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_abundance

FUNCTION_DESCRIPTION: This function calculates the effective Boltzmann constant and ion mass given some He/H abundance.

INPUTS:

OUTPUTS:

***********************************************************************************/

void ebtel_calc_abundance(void)
{
	K_B = 1.38e-16;
	double m_p = 1.67e-24;
	double m_p_old = m_p;

	//Calculate average ion mass
    double n_he_n_p = 0.075;   //He/p abundance.
    //Z_AVG = (1.0 + 2.0*n_he_n_p)/(1.0 + n_he_n_p); //Include Helium
    Z_AVG = 1.; //For Hydrad comparison.
    KB_FACT = 0.5*(1.0+1.0/Z_AVG); //modify equation of state for non-e-p plasma
    //double m_fact = (1.0 + n_he_n_p*4.0)/(2.0 + 3.0*n_he_n_p); //Include Helium
    double m_fact = (1 + n_he_n_p*4.)/2.; //For Hydrad comparison

	//Set global variables
    M_P = m_p*m_fact*(1.0 + Z_AVG)/Z_AVG; //Average ion mass
	MU = M_P/m_p_old;	//Mean molecular weight
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_ic

FUNCTION_DESCRIPTION: This function calculates the initial temperature using the static equilibrium
conditions for the hydrodynamic equations.

INPUTS:
	kpar--temperature bins used to calculate the radiative loss function
	r3--ratio of TR and coronal radiative losses (coefficient c1 from Papers I,II)
	loop_length--loop half-length (Mm)
	opt--Option structure of input values

OUTPUTS:
	return_array--array that holds r3, rad, t, n, p, and v

***********************************************************************************/

double * ebtel_calc_ic(double r3, double loop_length, struct Option *opt)
{
	//Variable declarations for both cases
	double *return_array = malloc(sizeof(double[6]));
	double r2 = ebtel_calc_c2();
	double heat = ebtel_heating(0,opt);

	if(strcmp(opt->ic_mode,"force") == 0 || strcmp(opt->ic_mode,"st_eq") == 0)
	{
		//Variable declarations
		int i;
		double tt_old;
		double nn_old;
		double tt_new;
		double rad;
		double nn;
		double p;
		double v;
		double err;
		double err_n;
		double tol;

		//Check if the heating array begins with a zero. If so, return an error.
		if (heat == 0.)
		{
			printf("ERROR! No initial loop heating: heat(0)=0. Provide valid heating input. Exiting the program\n");
			exit(0);
		}

		//First set up trial values for static equilibrium (i.e. d/dt = 0)
		tt_old = r2*pow(3.5*r3/(1 + r3)*loop_length*loop_length*heat/KAPPA_0_E,TWO_SEVENTHS);
		rad = ebtel_rad_loss(tt_old,opt->rad_option);
		nn = pow(heat/((1+r3)*rad),0.5);
		nn_old = nn;

		//Compute initial values for parameters t and n by iterating on temperature (tt) and R_tr/R_c (r3)
		tol = 1e+3;		//error tolerance

		for(i=0; i<=100; i++)
		{
			r3 = ebtel_calc_c1(tt_old,tt_old,nn,loop_length,rad,opt);										//recalculate r3 coefficient
			tt_new = r2*pow((3.5*r3/(1+r3)*pow(loop_length,2)*heat/KAPPA_0_E),TWO_SEVENTHS);	//temperature at new r3
			rad = ebtel_rad_loss(tt_new,opt->rad_option);											//radiative loss at new temperature
			nn = pow(heat/((1+r3)*rad),0.5);												//density at new r3 and new rad
			err = tt_new - tt_old;															//difference between t_i, T_i-1
			err_n = nn - nn_old;
			//Break the loop if the error gets below a certain threshold
			if(fabs(err)<tol && fabs(err_n)<tol)
			{
				tt_old = tt_new;
				nn_old = nn;
				break;
			}
			tt_old = tt_new;
			nn_old = nn;
		}

		//Calculate the density
		nn = pow(heat/((1+r3)*rad),0.5);

		//To use parameters consistent with the cases invoked in Paper II, we read in initial values for n,T rather than
		//calculating them using scaling laws or static equilibrium
		if(strcmp(opt->ic_mode,"force") == 0)
		{
			tt_old = opt->T0;
			nn = opt->n0;

			//DEBUG
			//These quantities being read in are actually apex quantities so we need to account for this
			//If later on forced ICs are not apex quantities, do not use this portion of the code.
			/*
			tt_old = r2*tt_old;
			double sc = ebtel_calc_lambda(tt_old);
			nn = nn/r2*exp(2*loop_length/(PI*sc)*(1. - sin(PI/5.)));
			*/
		}

		//Calculate resulting pressure, velocity
		p = KB_FACT*K_B*nn*tt_old;
		v = 0.;

		//Set array values
		return_array[0] = r3;
		return_array[1] = rad;
		return_array[2] = tt_old;
		return_array[3] = nn;
		return_array[4] = p;
		return_array[5] = v;
	}
	else if(strcmp(opt->ic_mode,"scaling") == 0)
	{
		//Variable declarations
		double lambda_0;
		double bb;
		double t_0;
		double p_0;
		double n_0;
		double v_0;
		double rad;

		//Alternatively, we could use the scaling laws to determine our initial conditions
		lambda_0 = 1.95e-18;			//lambda = lambda_0*T
		bb = -TWO_THIRDS;//-0.5			//power law for radiative loss function
		t_0 = r2*pow((3.5/KAPPA_0_E*heat),TWO_SEVENTHS)*pow(loop_length,2.0*TWO_SEVENTHS);
		p_0 = pow(r2,-SEVEN_HALVES*0.5)*pow(8.0/7.0*KAPPA_0_E/lambda_0,0.5)*KB_FACT*K_B*pow(t_0,((11.0-2.0*bb)/4.0))/loop_length;
		n_0 = 0.5*p_0/(KB_FACT*K_B*t_0);
		v_0 = 0;

		//Set array values
		rad = ebtel_rad_loss(t_0,opt->rad_option);
		return_array[0] = ebtel_calc_c1(t_0,t_0,n_0,loop_length,rad,opt);
		return_array[1] = rad;
		return_array[2] = t_0;
		return_array[3] = n_0;
		return_array[4] = p_0;
		return_array[5] = v_0;
	}
	else
	{
		printf("Invalid heat flux calculation option.\n");
		exit(0);
	}


	return return_array;
}


/***********************************************************************************

FUNCTION NAME: ebtel_calc_conduction

FUNCTION_DESCRIPTION: This function calculates the thermal conduction for some temperature
T. Depending on whether a classical or dynamic calculation is selected, the function also
includes flux limiting.

INPUTS:
	T_e--electron temperature (K)
	T_i--ion temperature (K)
	n--number density (cm^-3)
	L--loop half-length (cm)
	flux_key--tell whether to use classical (0) or dynamic (1) calculation

OUTPUTS:
	F_e--electron heat flux
	F_i--ion heat flux
	F_eq--equilibrium heat flux

***********************************************************************************/

double * ebtel_calc_conduction(double T_e, double T_i, double n, double L, double rad, double r3, double sat_limit, char *heat_flux_option)
{

	//Declare variables
	double c1_e, c1_i;
	double c_sat;
	double f_cl_e, f_cl_i;
	double f_sat_e, f_sat_i;
	double f_e, f_i;
	double f_eq;
	double r2 = ebtel_calc_c2();
	double *flux_ptr = malloc(sizeof(double[3]));

	//Set up thermal conduction parameters (NEED TO CHANGE FOR e- AND ion)
	c1_e = -TWO_SEVENTHS*KAPPA_0_E;
	c1_i = -TWO_SEVENTHS*KAPPA_0_I;
	c_sat = -1.5*pow(K_B,1.5)/pow(M_EL,0.5);

	//Set up thermal conduction at the base
	f_cl_e = c1_e*pow(T_e/r2,SEVEN_HALVES)/L;	//Classical heat flux calculation
	f_cl_i = c1_i*pow(T_i/r2,SEVEN_HALVES)/L;

	//Decide on whether to use classical or dynamic heat flux
	if(strcmp(heat_flux_option,"classical")==0)
	{
		f_e = f_cl_e;
		f_i = f_cl_i;
	}
	else if(strcmp(heat_flux_option,"limited")==0)
	{
		//Compute flux limit
		f_sat_e = sat_limit*c_sat*n*pow(T_e,1.5);
		f_sat_i = sat_limit*c_sat*n*pow(T_i,1.5);

		//Compute final flux value
		f_e = -f_cl_e*f_sat_e/pow((pow(f_cl_e,2.) + pow(f_sat_e,2)),0.5);
		f_i = -f_cl_i*f_sat_i/pow((pow(f_cl_i,2.) + pow(f_sat_i,2)),0.5);

	}
	else
	{
		printf("Invalid heat flux option.\n");
		exit(0);
	}

	//Calculate equilibrium thermal conduction at base (-R_tr in Paper I)
	f_eq = -r3*pow(n,2.)*rad*L;

	//Set the flux array
	flux_ptr[0] = f_e;
	flux_ptr[1] = f_i;
	flux_ptr[2] = f_eq;

	//Return the array
	return flux_ptr;

}


/***********************************************************************************

FUNCTION NAME: ebtel_collision_freq

FUNCTION_DESCRIPTION: This function calculates the coulomb collision frequency between
the ions and electrons as given in Braginskii (1965) and used in Bradshaw and Cargill (2013). The
coulomb logarithm is calculated using the procedure found in Physics of the Solar Corona by
M.J. Aschwanden.

INPUTS:
	T_e--electron temperature (K)
	T_i--ion temperature (K)
	n--number density (cm^-3)

OUTPUTS:
	nu_ei--electron-ion collision frequency

***********************************************************************************/

double ebtel_collision_freq(double T_e, double T_i, double n)
{
	//Declare variables
	double ln_lambda;
	double nu_ei;
	double beta_1 = 1.0e+13;
	double beta_2 = 1.602*1e-9;

	//Expression for the Coulomb logarithm from Physics of the Solar Corona by M.J. Aschwanden
	ln_lambda = 23. - log(sqrt(n/beta_1)*pow(K_B*T_e/beta_2,-3./2.));

	//Calculate collision frequency
	nu_ei = 16./3.*sqrt(PI)*pow(Q_E,4.)/(M_EL*M_P)*pow(2.*K_B*T_e/M_EL,-3./2.)*n*ln_lambda;

	return nu_ei;
}


/***********************************************************************************

FUNCTION NAME: ebtel_calc_vel

FUNCTION_DESCRIPTION: This function calculates flow velocity. It is used outside of
the solver to update the velocity.

INPUTS:
	T_e--electron temperature
	T_i--ion temperature
	f_e--electron conductive flux
	f_i--ion conductive flux
	f_eq--equilibrium conductive flux
	par--parameter structure

OUTPUTS:
	v--flow velocity (cm/s)

***********************************************************************************/

double ebtel_calc_vel(double T_e, double T_i, double p_e, struct rk_params par)
{
	//Declare variables
	double r2e,r2i,r1e,r1i;
	double xi;
	double R_tr;
	double psi;
	double p_ev;
	double v;

	//Calculate ratio of base temperatures for ions and electrons
	r2e = ebtel_calc_c2();
	r2i = ebtel_calc_c2();
	r1e = ebtel_calc_c3();
	r1i = ebtel_calc_c3();
	xi = r1e/r1i*r2i/r2e*T_e/T_i;

	//Calculate the radiative loss of the transition region
	R_tr = -par.f_eq;

	//Approximate TR integral of v*dPe/ds term
	psi = (par.f_e - xi*par.f_i + R_tr)/(1 + xi);

	//Calculate enthalpy flux
	p_ev = 2./5.*(psi - par.f_e - R_tr);

	//Calculate v
	v = p_ev/p_e*par.r4;

	return v;
}

/***********************************************************************************

FUNCTION NAME: ebtel_thermal_conduction_timescale

FUNCTION_DESCRIPTION: This function calculates the minimum thermal conduction timescale
according to the hydrodynamic energy equation.

INPUTS:
	T_e--electron temperature
	T_i--ion temperature
	n--number density
	opt--structure containing input options

OUTPUTS:
	tau_tc--thermal conduction timescale

***********************************************************************************/

double ebtel_thermal_conduction_timescale(double T_e, double T_i, double n, struct Option *opt)
{
	//Variable declarations
	double temp;
	double kappa;

	//First determine whether T_e or T_i is higher
	if(T_e > T_i)
	{
		temp = T_e;
		kappa = KAPPA_0_E;
	}
	else
	{
		temp = T_i;
		kappa = KAPPA_0_I;
	}

	return 3.0*K_B*n*pow(((opt->loop_length)*1.0e+8),2.0)/(kappa*pow(temp,2.5));
}
