/***********************************************************************************

FILENAME: ebtel_functions_util.c

AUTHOR Will Barnes

DATE: created: 7 April 2014

DESCRIPTION: This file contains utility functions needed in EBTEL. For additional 
details, see ebtel_main.c.

***********************************************************************************/

//Include appropriate header file
#include "ebtel-2fl_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_print_header

FUNCTION_DESCRIPTION: This function prints a header to the screen and tells the user
what input is being used in this implementation of the model. 

INPUTS:
	n--total number of steps
	heating_shape--type of heating function used
	loop_length--half-length of the coronal loop (Mm)
	total_time--time the simulation was run for (s)
	Option opt--data structure that holds input parameters

OUTPUTS:
	

***********************************************************************************/

void ebtel_print_header(int n, struct Option *opt)
{
	//Print a header and tell the user what options are being used to begin the model
	printf("************************************************************************************\n");
	printf("            Enthalpy Based Thermal Evolution of Loops (EBTEL)						\n");
	printf("			   Two-fluid Model											\n");
	printf("************************************************************************************\n\n");
	printf("Original single-fluid code written in IDL by J.A. Klimchuk, S. Patsourakos, P.J. Cargill\n");
	printf("See Klimchuk, J.A, S. Patsourakos & P.J. Cargill 2008, ApJ 682:1351-2362\n");
	printf("See also Cargill, P.J., S.J. Bradshaw & J.A. Klimchuk 2012, ApJ 752:161-174\n\n");
	printf("Translation into the C Programming Language by Will Barnes,\nDept. of Physics & Astronomy, Rice University (2014)\n");
	printf("************************************************************************************\n\n");
	printf("INPUTS\n");
	printf("------\n");
	printf("Total time: %d s\n",opt->total_time);
	printf("Time step: %f s\n",opt->tau);
	printf("Loop half-length: %f Mm\n",opt->loop_length);
	printf("Usage option(see documentation): %s\n",opt->usage_option);
	printf("Heating pulse shape: %s\n",opt->heating_shape);
	printf("Heating species: %s\n",opt->heat_species);
	if(strcmp(opt->solver,"rk4")==0)
	{printf("Solving equations using fourth order Runge-Kutta routine\n");
	}
	else if(strcmp(opt->solver,"rka4")==0)
	{printf("Solving equations using adaptive fourth order Runge-Kutta routine\n");
	}
	else if(strcmp(opt->solver,"euler")==0)
	{printf("Solving equations using Euler method\n");
	}
	else
	{printf("Invalid solver option\n");
	}
	if(strcmp(opt->rad_option,"rtv")==0)
	{printf("Using Rosner-Tucker-Vaiana Loss Function\n");
	}
	else if(strcmp(opt->rad_option,"rk")==0)
	{printf("Using Raymond-Klimchuk Loss Function\n");
	}
	else
	{printf("Invalid radiative loss option\n");
	}
	printf("Using %s method to calculate the heat flux\n",opt->heat_flux_option);
	if(strcmp(opt->usage_option,"dem")==0 || strcmp(opt->usage_option,"rad_ratio")==0)
	{
		printf("Using %s method to calculate DEM in the TR\n",opt->dem_option);
	}
	if(strcmp(opt->ic_mode,"st_eq")==0)
	{printf("Using static equilibrium to calculate initial conditions\n");
	}
	else if(strcmp(opt->ic_mode,"force")==0)
	{printf("Forcing initial conditions with T_0 = %f MK and n_0 = %f*10^8 cm^-3\n",opt->T0/pow(10,6),opt->n0/pow(10,8));
	}
	else if(strcmp(opt->ic_mode,"scaling")==0)
	{printf("Using scaling laws to calculate initial conditions\n");
	}
	else
	{printf("Invalid initial conditions option\n");
	}
	printf("\n");
}

/***********************************************************************************

FUNCTION NAME: ebtel_input_setter

FUNCTION_DESCRIPTION: This function reads the relevant inputs from an XML config file 
and sets them to a structure. The xmllib2 library is used. The necessary compiler flags 
are set in the makefile or can be configured manually.


INPUTS:
		char filename--name of xml file that holds input configuration
OUTPUTS:
		struct Option opt--pointer to structure that holds input configuration
	
***********************************************************************************/

struct Option *ebtel_input_setter(char *filename)
{
	//Declare temp char pointer to hold results
	char *temp;
	
	//Declare doc and root pointers
	xmlDocPtr doc;
	xmlNodePtr root;
	//Create the document tree
	doc = xmlParseFile(filename); 
	//Check if the document can be parsed
	if(doc == NULL)
	{
		printf("%s cannot be parsed or does not exist. Please specify another file.\n",filename);
		exit(0);
	}
	//Point root at first child of the tree
	root = doc->children;
	
	//Declare the structure that will be returned
	struct Option *opt = malloc(sizeof(struct Option));
	assert(opt != NULL);
	
	//Set the value of each structure field from the xml file
	//Double
	temp = ebtel_xml_reader(root,"total_time",NULL);
	opt->total_time = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"tau",NULL);
	opt->tau = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"loop_length",NULL);
	opt->loop_length = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"sat_limit",NULL);
	opt->sat_limit = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"h_nano",NULL);
	opt->h_nano = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"t_pulse_half",NULL);
	opt->t_pulse_half = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"t_start",NULL);
	opt->t_start = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"rka_error",NULL);
	opt->rka_error = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"T0",NULL);
	opt->T0 = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"n0",NULL);
	opt->n0 = atof(temp);
	free(temp);
	temp = NULL;

	temp = ebtel_xml_reader(root,"h_back",NULL);
	opt->h_back = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"mean_t_start",NULL);
	opt->mean_t_start = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"std_t_start",NULL);
	opt->std_t_start = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"alpha",NULL);
	opt->alpha = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"amp0",NULL);
	opt->amp0 = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"amp1",NULL);
	opt->amp1 = atof(temp);
	free(temp);
	temp = NULL;
	
	//Int
	temp = ebtel_xml_reader(root,"index_dem",NULL);
	opt->index_dem = atoi(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"num_events",NULL);
	opt->num_events = atoi(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"sample_rate",NULL);
	opt->sample_rate = atoi(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"r3_rad_0",NULL);
	opt->r3_rad_0 = atof(temp);
	free(temp);
	temp = NULL;
	
	temp = ebtel_xml_reader(root,"r3_cond_0",NULL);
	opt->r3_cond_0 = atof(temp);
	free(temp);
	temp = NULL;
	
	//Char
	opt->heating_shape = ebtel_xml_reader(root,"heating_shape",NULL);
	opt->usage_option = ebtel_xml_reader(root,"usage_option",NULL);
	opt->rad_option = ebtel_xml_reader(root,"rad_option",NULL);
	opt->dem_option = ebtel_xml_reader(root,"dem_option",NULL);
	opt->heat_flux_option = ebtel_xml_reader(root,"heat_flux_option",NULL);
	opt->solver = ebtel_xml_reader(root,"solver",NULL);
	opt->ic_mode = ebtel_xml_reader(root,"ic_mode",NULL);
	opt->output_file = ebtel_xml_reader(root,"output_file",NULL);
	opt->print_plasma_params = ebtel_xml_reader(root,"print_plasma_params",NULL);
	opt->heat_species = ebtel_xml_reader(root,"heat_species",NULL);
	opt->t_start_switch = ebtel_xml_reader(root,"t_start_switch",NULL);
	opt->amp_switch = ebtel_xml_reader(root,"amp_switch",NULL);
	opt->t_end_switch = ebtel_xml_reader(root,"t_end_switch",NULL);
	opt->r3_loss_correction = ebtel_xml_reader(root,"r3_loss_correction",NULL);
	opt->r3_grav_correction = ebtel_xml_reader(root,"r3_grav_correction",NULL);
	
	//Free the document tree
	xmlFreeDoc(doc);
	
	//Return the structure
	return opt;	
}

/***********************************************************************************

FUNCTION NAME: ebtel_xml_reader

FUNCTION_DESCRIPTION: This function reads the relevant inputs from an XML config file 
recursively and is capable of processing XML files of arbitrary depth. When calling the
function, nodeValue should be initialized with NULL and root with the doc->children where
doc is the XML file document tree.

INPUTS:
		xmlNodePtr root -- root node to start the read
		char *nodeName -- name of the node we are searching for
		char *nodeValue -- value to retrieve (default: NULL)
OUTPUTS:
		char *nodeValue -- value to retrieve
	
***********************************************************************************/

char *ebtel_xml_reader(xmlNodePtr root, char *nodeName, char *nodeValue)
{
	//Declare new node instance
	xmlNodePtr cur;
	//Point it to the child of the root passed to the function
	cur = root->children;
	
	//Begin loop
	//Only break when the list is over or we have found the value
	while(cur != NULL && nodeValue == NULL)
	{
		//Check if the current node has children
		//If not, cur is a value so we check it
		if(cur->children == NULL)
		{
			//Check the parent of the node
			if(!xmlIsBlankNode(cur) && strcmp((char*)root->name,nodeName)==0)
			{
				//Get the value of the node
				nodeValue = (char *)xmlNodeGetContent(cur);
				//Return the value
				return nodeValue;
			}
		}
		
		//If the current node has children, descend a level deeper to check values 
		else
		{
			//Call recursively
			nodeValue = ebtel_xml_reader(cur,nodeName,nodeValue);
			//Break the loop if we have found the value
			if(nodeValue != NULL)
			{
				break;
			}
		}
		//Move to the next node
		cur = cur->next;
	}
	
	//Return the value
	return nodeValue;
}

/***********************************************************************************

FUNCTION NAME: ebtel_file_writer

FUNCTION_DESCRIPTION: This function prints the data from the structure ebtel_params_st
to a file. It produces separate files for the loop data and DEM data. 

INPUTS:
	loop_length--loop half length in Mm
	n--length of the arrays to be printed
	opt--Option structure with inputs necessary for labeling
	params_final--pointer to structure of output data
	
OUTPUTS:
	none

***********************************************************************************/

void ebtel_file_writer(struct Option *opt, struct ebtel_params_st *params_final)
{
	//Declare variables
	int i;
	int n = params_final->i_max;
	char filename_out[250];
	char filename_out_dem[250];
	char filename_out_heat_amp[250];
	FILE *out_file;
	
	if(strcmp(opt->print_plasma_params,"True") == 0)
	{
		//Use either the default filename or a custom name specified in the configuration file
		if(strcmp(opt->output_file,"default") == 0)
		{
			//Check and see if directory 'data' exists. If it does not, then create a new one.
			struct stat st = {0};
			if(stat("../data",&st) == -1){
				mkdir("../data",0777);
			}
		
			//Open the file that we are going to write the data to 
			sprintf(filename_out,"../data/ebteldatL%.*f_%s_%s_%s.txt",1,opt->loop_length,opt->usage_option,opt->heating_shape,opt->solver);	
			
		}
		else
		{
			sprintf(filename_out,"%s.txt",opt->output_file);
		}
		
		//Open the output stream
		out_file = fopen(filename_out,"wt");
		
		//Check to make sure the file was opened successfully
		if(out_file == NULL)
		{
			printf("The file %s could not be opened for writing.\n",filename_out);
			exit(0);
		}
		
		//The members of the structure params_final have now been set so we need to unpack them and set our arrays so that we can easily save our data.
		for(i = 0; i<n; i++)
		{	
			//Only save every sample_rate timestep 
			if(i % opt->sample_rate == 0)
			{
				//If we used usage_option 'rad_ratio', then we need to save f_ratio and rad_ratio as well.
				//We set them to zero otherwise just as a placeholder.
				if(strcmp(opt->usage_option,"rad_ratio") == 0)
				{	
					//Print the data to the file filename using tab delimited entries
					fprintf(out_file,"%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",*(params_final->time + i),*(params_final->temp_e + i),*(params_final->temp_i + i),*(params_final->ndens + i),*(params_final->press_e + i),*(params_final->press_i + i),*(params_final->vel + i),*(params_final->tapex_e + i),*(params_final->tapex_i + i),*(params_final->napex +i),*(params_final->papex_e + i),*(params_final->papex_i + i),*(params_final->cond_e + i),*(params_final->cond_i + i),*(params_final->rad_cor + i),*(params_final->heat + i),*(params_final->coeff_1 + i),*(params_final->rad + i),*(params_final->tau + i),*(params_final->rad_ratio + i),*(params_final->f_ratio + i));
				}
				else
				{
					//Print the data to the file filename using tab delimited entries
					fprintf(out_file,"%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",*(params_final->time + i),*(params_final->temp_e + i),*(params_final->temp_i + i),*(params_final->ndens + i),*(params_final->press_e + i),*(params_final->press_i + i),*(params_final->vel + i),*(params_final->tapex_e + i),*(params_final->tapex_i + i),*(params_final->napex +i),*(params_final->papex_e + i),*(params_final->papex_i + i),*(params_final->cond_e + i),*(params_final->cond_i + i),*(params_final->rad_cor + i),*(params_final->heat + i),*(params_final->coeff_1 + i),*(params_final->rad + i),*(params_final->tau + i));
				}
			}
			
		}
		
		//Close the file when we are done writing data to it
		fclose(out_file);
		
		//Tell the user where the results were printed
		printf("The results were printed to the file %s\n",filename_out);	
	}

	
	//If we used the power-law or file option, print the resulting amplitudes to a file
	if(strcmp(opt->amp_switch,"power_law")==0 || strcmp(opt->amp_switch,"file")==0)
	{
		if(strcmp(opt->output_file,"default") == 0)
		{
			sprintf(filename_out_heat_amp,"../data/ebtel-2fldatL%.*f_%s_%s_%s_heat_amp.txt",1,opt->loop_length,opt->usage_option,opt->heating_shape,opt->solver);
		}
		else
		{
			sprintf(filename_out_heat_amp,"%s_heat_amp.txt",opt->output_file);
		}
		
		//Open the heating amplitude file
		out_file = fopen(filename_out_heat_amp,"wt");
		
		//Check to make sure the file was opened successfully
		if(out_file == NULL)
		{
			printf("The file %s could not be opened for writing.\n",filename_out_heat_amp);
			exit(0);
		}
		
		for(i = 0; i<opt->num_events; i++)
		{
			//Write the data to a file
			fprintf(out_file,"%f\n",*(opt->amp + i));
		}
		
		//Close the heating amplitude file
		fclose(out_file);
		
		//Tell the user where the results were printed
		printf("The heating amplitude results were printed to the file %s\n",filename_out_heat_amp);
	}
	
	//If we chose to calculate the TR DEM, we need to write this data to a separate file.
	if(strcmp(opt->usage_option,"dem")==0 || strcmp(opt->usage_option,"rad_ratio")==0)
	{
		//Make the DEM data filename--check for a custom output filename
		if(strcmp(opt->output_file,"default")==0)
		{
			sprintf(filename_out_dem,"../data/ebtel-2fldatL%.*f_%s_%s_%s_dem.txt",1,opt->loop_length,opt->usage_option,opt->heating_shape,opt->solver);			
		}
		else
		{
			sprintf(filename_out_dem,"%s_dem.txt",opt->output_file);
		}
		
		//Open the DEM data file
		out_file = fopen(filename_out_dem,"wt");
		
		//Check to make sure the file was opened successfully
		if(out_file == NULL)
		{
			printf("The file %s could not be opened for writing.\n",filename_out_dem);
			exit(0);
		}
		
		for(i=0; i<opt->index_dem; i++)
		{	
			//Now write this data to a file. 
			fprintf(out_file,"%e\t%e\t%e\t%e\t%e\n",*(params_final->logtdem + i), *(params_final->dem_tr_log10mean + i),  *(params_final->dem_cor_log10mean + i), *(params_final->dem_tot_log10mean + i), *(params_final->em_cor_log10mean + i));
		}
		
		//Close the DEM data file
		fclose(out_file);
		
		//Tell the user where the DEM data was printed to
		printf("The DEM results were printed to the file %s\n",filename_out_dem);
	}
	
}

/***********************************************************************************

FUNCTION NAME: ebtel_linspace

FUNCTION_DESCRIPTION: This function is similar to the linspace function in MATLAB. It creates a one-dimensional array with specified number of entries and starting and ending point. 

INPUTS:
	a--starting point
	b--ending point
	n--number of entries
	
OUTPUTS:
	linspace--pointer to array with n points starting at a and ending at b

***********************************************************************************/

double * ebtel_linspace( int a, int b, int n)
{
	//Declare necessary variables
	int i;
	double B = b;	//cast as doubles
	double A = a;	
	double N = n;	
	double *linspace = malloc(sizeof(double[n]));
	double interval = (B - A)/(N-1);	//spacing between points
	
	//Make the array
	linspace[0] = a;
	for(i = 1; i<n; i++)
	{
		linspace[i] = interval + *(linspace + (i-1));
	}

	return linspace;
}

/***********************************************************************************

FUNCTION NAME: ebtel_colon_operator

FUNCTION_DESCRIPTION: This function is similar to using the colon operator in MATLAB. It takes the starting point, end point, and spacing and creates an array beginning at a and ending at b with each of the points separated by d. If it cannot make it exactly to b, the function will always undershoot b. 

INPUTS:
	a--starting point
	b--ending point
	d--spacing between entries
	
OUTPUTS:
	colon_array--pointer to array with values between a and b spaced by d

***********************************************************************************/

double * ebtel_colon_operator(double a, double b, double d)
{
	//Declare variables
	int i = 0;					//initialize counter
	int size = ceil((b - a)/d);	//set array size
	if(size == (b-a)/d)			//add 1 to size if there is no underflow in the array
	{
		size = size + 1;
	}
	double *colon_operator = malloc(sizeof(double[size]));	//reserve memory for the pointer
	
	//Begin the loop to set the array
	while(*(colon_operator + i) < b)
	{
		colon_operator[i] = a + i*d;
		i++;
		
		//If we've reached the condition, we don't want to set the next entry
		if(a + i*d > b)
		{
			break;
		}
	}
	
	//Return the pointer
	return colon_operator;
}


 /**********************************************************************************
 
 Function name: ebtel_avg_val
 
 Function description: This function accepts an array and calculates its mean or average
 value. 
 
 Input
	double numbers []	Array of double values
	int length			Integer value, length of array numbers
 
 Return
	double mean			Double value returned by the function.
 
 *********************************************************************************/
 
 double ebtel_avg_val (double numbers[], int length)
 {
 	//Declare some variables
 	int i;
 	double mean;
 	double sum = 0;
 	
 	//Calculate the sum
 	for (i = 0; i<length; i++)
 	{
 		sum = sum + numbers[i];	
 	}
 	
 	//Calculate the mean
 	mean = sum/length;
 	
 	//Return the mean
 	return mean;
 }
 
 /**********************************************************************************
 
 Function name: ebtel_weighted_avg_val
 
 Function description: This function accepts an array and calculates its weighted 
 average. 
 
 Input
	double numbers 		Array of double values
	int length			Integer value, length of array numbers
 	double weight		Array of weights for each entry in array numbers 
 						(Length must be the same as numbers)
 
 Return
	double mean			Double value returned by the function.
 
 *********************************************************************************/
 
 double ebtel_weighted_avg_val(double *numbers, int length, double *weight)
 {
  	//Declare some variables
  	int i;
  	double mean = 0.;
	double rel_weight;
  	double sum = 0.;
 	
  	//Calculate the sum of the weights
  	for (i = 0; i<length; i++)
  	{
  		sum = sum + *(weight+ i);	
  	}
	
	//Construct the average by calculating the respective weights and multiplying
	//each entry in numbers[] and then summing
	for(i=0;i<length;i++)
	{
		rel_weight = *(weight + i)/sum;
		mean = mean + *(numbers+ i)*rel_weight;
	}
	
	//Return the weighted mean
	return mean;
 }
 
 
/**********************************************************************************
 
 Function name: ebtel_max_val
 
 Function description: This function accepts two numbers and returns the maximum of 
 those two numbers.
 
 Input
	double num_1		first number to compare
	double num_2		second number to compare	
 
 Return
	double max_val		maximum of num_1 and num_2
 
 *********************************************************************************/
 
 double ebtel_max_val(double num_1, double num_2)
 {	
 	if(num_1 > num_2)	
 	{
 		return num_1;
 	}
 	else if(num_2 > num_1)
 	{
 		return num_2;
 	}
 	else  
 	{
 		return num_2;
 	}
 }
 
 /**********************************************************************************
 
 Function name: ebtel_min_val
 
 Function description: This function accepts two numbers and returns the minimum of
 these two numbers.
 
 Input
	double num_1		first number to compare
	double num_2		second number to compare	
 
 Return
	double min_val		minimum of num_1 and num_2
 
 *********************************************************************************/
 
 double ebtel_min_val(double num_1, double num_2)
 {
 	
 	if(num_1 < num_2)
 	{
 		return num_1;
 	}
 	else if(num_2 < num_1)
 	{
 		return num_2;
 	}
 	else
	{					
 		return num_1;
 	}
 }
 
 /**********************************************************************************
 
 FUNCTION NAME: ebtel_reallocate_mem
 
 FUNCTION DESCRIPTION: This function reallocates memory to the arrays in the structure
 struct ebtel_params_st as well as the two-dimensional arrays that store the coronal and
 TR DEM values.
 
 INPUT:
 	mem_lim			integer value of current memory limit
 	new_mem_lim		integer value of expanded memory limit
	par_struct		pointer to structure memory that will be freed

 OUTPUT:
	none
 
 *********************************************************************************/
 
 void ebtel_reallocate_mem(int mem_lim, int new_mem_lim, struct ebtel_params_st *par_struct, struct Option *opt)
 {
	 //Reallocate each of the arrays in par_struct
	 double *time_r;
	 if((time_r = realloc(par_struct->time,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->time = time_r;
	 
	 double *tau_r;
	 if((tau_r = realloc(par_struct->tau,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->tau = tau_r;
	 
	 double *heat_r;
	 if((heat_r = realloc(par_struct->heat,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->heat = heat_r;
	 
	 double *temp_e_r;
	 if((temp_e_r = realloc(par_struct->temp_e,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->temp_e = temp_e_r;
	 
	 double *temp_i_r;
	 if((temp_i_r = realloc(par_struct->temp_i,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->temp_i = temp_i_r;
	 
	 double *ndens_r;
	 if((ndens_r = realloc(par_struct->ndens,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->ndens = ndens_r;
	 
	 double *press_e_r;
	 if((press_e_r = realloc(par_struct->press_e,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->press_e = press_e_r;
	 
	 double *press_i_r;
	 if((press_i_r = realloc(par_struct->press_i,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->press_i = press_i_r;
	 
	 double *vel_r;
	 if((vel_r = realloc(par_struct->vel,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->vel = vel_r;
	 
	 double *tapex_e_r;
	 if((tapex_e_r = realloc(par_struct->tapex_e,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->tapex_e = tapex_e_r;
	 
	 double *tapex_i_r;
	 if((tapex_i_r = realloc(par_struct->tapex_i,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->tapex_i = tapex_i_r;
	 
	 double *napex_r;
	 if((napex_r = realloc(par_struct->napex,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->napex = napex_r;
	 
	 double *papex_e_r;
	 if((papex_e_r = realloc(par_struct->papex_e,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->papex_e = papex_e_r;
	 
	 double *papex_i_r;
	 if((papex_i_r = realloc(par_struct->papex_i,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->papex_i = papex_i_r;
	 
	 double *coeff_1_r;
	 if((coeff_1_r = realloc(par_struct->coeff_1,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->coeff_1 = coeff_1_r;
	 
	 double *cond_e_r;
	 if((cond_e_r = realloc(par_struct->cond_e,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->cond_e = cond_e_r;
	 
	 double *cond_i_r;
	 if((cond_i_r = realloc(par_struct->cond_i,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->cond_i = cond_i_r;
	 
	 double *rad_cor_r;
	 if((rad_cor_r = realloc(par_struct->rad_cor,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->rad_cor = rad_cor_r;
	 
	 double *rad_r;
	 if((rad_r = realloc(par_struct->rad,sizeof(double[new_mem_lim]))) == NULL)
	 {
		 printf("Reallocation failed. Out of memory.\n");
		 exit(0);
	 }
	 par_struct->rad = rad_r;
	 
	 //Check if rad_ratio or f_ratio were set
	 if(strcmp(opt->usage_option,"rad_ratio")==0)
	 {
		 double *f_ratio_r;
		 if((f_ratio_r = realloc(par_struct->f_ratio,sizeof(double[new_mem_lim]))) == NULL)
		 {
			 printf("Reallocation failed. Out of memory.\n");
			 exit(0);
		 }
		 par_struct->f_ratio = f_ratio_r;
		 
		 double *rad_ratio_r;
		 if((rad_ratio_r = realloc(par_struct->rad_ratio,sizeof(double[new_mem_lim]))) == NULL)
		 {
			 printf("Reallocation failed. Out of memory.\n");
			 exit(0);
		 }
		 par_struct->rad_ratio = rad_ratio_r;
	 }
	 
 }
 
 /**********************************************************************************
 
 FUNCTION NAME: ebtel_reallocate_two_d_array
 
 FUNCTION DESCRIPTION: This function reallocates memory for a doubly referenced pointer (i.e. a two-dimensional array)
 
 INPUT:
 	array_2d -- doubly referenced pointer
 	mem_lim -- current memory limit for the index that is being reallocated on
 	new_mem_lim -- new memory limit for the index that is being reallocated on
 	dim_2 -- static index
 
 OUTPUT:
 	array_2d -- reallocated doubly referenced pointer 
 
 *********************************************************************************/
 
 double **ebtel_reallocate_two_d_array(double **array_2d, int mem_lim, int new_mem_lim, int dim_2)
 {
	int k;
	double **tmp;
	if((tmp = (double **)realloc(array_2d,sizeof(double *)*new_mem_lim)) == NULL )
	{
		printf("Out of memory!\n");
		exit(0);
	}
	for(k = mem_lim; k<new_mem_lim; k++)
	{
		tmp[k] = (double *)malloc(sizeof(double)*dim_2);
	}

	array_2d = tmp;

	return array_2d;
 } 

 /**********************************************************************************
 
 Function name: ebtel_free_mem
 
 Function description: This function frees up memory that was malloc'd by the ebtel_loop
 _solver function in order to return the param structure
 
 Input
	struct ebtel_params_st *par_struct		pointer to structure memory that will be freed	
 
 Return
	none
 
 *********************************************************************************/
 void ebtel_free_mem(struct ebtel_params_st *par_struct, struct Option *opt)
 {
 	//First check that the pointer is valid
 	assert(par_struct != NULL);
 	
 	//Free the memory of each of the structure members as well as the structure itself
 	free(par_struct->time);
 	par_struct->time = NULL;
	free(par_struct->tau);
	par_struct->tau = NULL;
 	free(par_struct->heat);
 	par_struct->heat = NULL;
	free(par_struct->coeff_1);
	par_struct->coeff_1 = NULL;
	
 	free(par_struct->temp_e);
 	par_struct->temp_e = NULL;
 	free(par_struct->temp_i);
 	par_struct->temp_i = NULL;
	free(par_struct->tapex_e);
	par_struct->tapex_e = NULL;
	free(par_struct->tapex_i);
	par_struct->tapex_i = NULL;
	
	free(par_struct->press_e);
	par_struct->press_e = NULL;
	free(par_struct->press_i);
	par_struct->press_i = NULL;
	free(par_struct->papex_e);
	par_struct->papex_e = NULL;
	free(par_struct->papex_i);
	par_struct->papex_i = NULL;
	
	free(par_struct->ndens);
	par_struct->ndens = NULL;
	free(par_struct->napex);
	par_struct->napex = NULL;
	free(par_struct->vel);
	par_struct->vel = NULL;
	
	free(par_struct->cond_e);
	par_struct->cond_e = NULL;
	free(par_struct->cond_i);
	par_struct->cond_i = NULL;
	free(par_struct->rad_cor);
	par_struct->rad_cor = NULL;
	free(par_struct->rad);
	par_struct->rad = NULL;
	
	//Free memory based on usage option
	if(strcmp(opt->usage_option,"dem")==0 || strcmp(opt->usage_option,"rad_ratio")==0)
	{
		if(strcmp(opt->usage_option,"rad_ratio")==0)
		{
			free(par_struct->f_ratio);
			par_struct->f_ratio = NULL;
			free(par_struct->rad_ratio);
			par_struct->rad_ratio = NULL;
		}
		
		free(par_struct->logtdem);
		par_struct->logtdem = NULL;
		free(par_struct->dem_tr_log10mean);
		par_struct->dem_tr_log10mean = NULL;
		free(par_struct->dem_cor_log10mean);
		par_struct->dem_cor_log10mean = NULL;
		free(par_struct->em_cor_log10mean);
		par_struct->em_cor_log10mean = NULL;
		free(par_struct->dem_tot_log10mean);
		par_struct->dem_tot_log10mean = NULL;
	}
	
	//Free memory reserved for the par structure
	free(par_struct);
	
	//Free memory used by char array members of opt
	free(opt->heat_flux_option);
	opt->heat_flux_option = NULL;
	free(opt->dem_option);
	opt->dem_option = NULL;
	free(opt->rad_option);
	opt->rad_option = NULL;
	free(opt->usage_option);
	opt->usage_option = NULL;
	free(opt->solver);
	opt->solver = NULL;
	free(opt->ic_mode);
	opt->ic_mode = NULL;
	free(opt->output_file);
	opt->output_file = NULL;
	free(opt->print_plasma_params);
	opt->print_plasma_params = NULL;
	free(opt->heat_species);
	opt->heat_species = NULL;	
	free(opt->heating_shape);
	opt->heating_shape = NULL;
	free(opt->t_start_switch);
	opt->t_start_switch = NULL;
	free(opt->t_end_switch);
	opt->t_end_switch = NULL;	
	free(opt->amp_switch);
	opt->amp_switch = NULL;
	free(opt->r3_loss_correction);
	opt->r3_loss_correction = NULL;
	free(opt->r3_grav_correction);
	opt->r3_grav_correction = NULL;
	
	//Free the t_start, amp, and t_end arrays
	free(opt->t_start_array);
	opt->t_start_array = NULL;
	free(opt->amp);
	opt->amp = NULL;
	free(opt->t_end_array);
	opt->t_end_array = NULL;
	
	//Free memory used by opt input structure
	free(opt);
 }
 