/*
loop.cpp
Loop object that will hold all information about the loop and be evolved in time.
*/

#include "loop.h"

Loop::Loop(char *ebtel_config, char *rad_config)
{
  tinyxml2::XMLDocument doc;
  tinyxml2::XMLElement *root;
  tinyxml2::XMLElement *heating;
  double helium_to_hydrogen_ratio;

  //Open file
  tinyxml2::XMLError load_ok = doc.LoadFile(ebtel_config);
  if(load_ok != 0)
  {
  	printf("Failed to load XML configuration file %s.\n",ebtel_config);
  	//TODO: Exit or break out from here
  }
  //Parse file and read into data structure
  root = doc.FirstChildElement();
  //Numeric parameters
  parameters.total_time = std::stod(get_element_text(root,"total_time"));
  parameters.tau = std::stod(get_element_text(root,"tau"));
  parameters.loop_length = std::stod(get_element_text(root,"loop_length"));
  parameters.rka_error = std::stod(get_element_text(root,"loop_length"));
  parameters.saturation_limit = std::stod(get_element_text(root,"saturation_limit"));
  parameters.c1_cond0 = std::stod(get_element_text(root,"c1_cond0"));
  parameters.c1_rad0 = std::stod(get_element_text(root,"c1_rad0"));
  helium_to_hydrogen_ratio = std::stod(get_element_text(root,"helium_to_hydrogen_ratio"));
  //Boolean parameters
  parameters.use_c1_loss_correction = string2bool(get_element_text(root,"use_c1_loss_correction"));
  parameters.use_c1_grav_correction = string2bool(get_element_text(root,"use_c1_grav_correction"));
  parameters.use_power_law_radiative_losses = string2bool(get_element_text(root,"use_power_law_radiative_losses"));
  parameters.use_spitzer_conductivity = string2bool(get_element_text(root,"use_spitzer_conductivity"));
  parameters.calculate_dem = string2bool(get_element_text(root,"calculate_dem"));
  //String parameters
  parameters.solver = get_element_text(root,"solver");
  parameters.output_filename = get_element_text(root,"output_filename");

  //Initialize radiation model object
  if(parameters.use_power_law_radiative_losses)
  {
    radiation_model = new CRadiation();
  }
  else
  {
    radiation_model = new CRadiation(rad_config,false);
  }

  //Initialize heating object
  heating = get_element(root,"heating");
  heater = new Heater(heating);

  //Initialize DEM object
  //TODO: implement after everything else is in place; only instantiate if we want to do that calculation

  doc.Clear();

  // Calculate needed He abundance corrections
  CalculateAbundanceCorrection(helium_to_hydrogen_ratio);

  //Estimate results array length
  N = int(ceil(parameters.total_time/parameters.tau));
  //Reserve memory for results
  results.time.resize(N);
  results.heat.resize(N);
  results.pressure_e.resize(N);
  results.pressure_i.resize(N);
  results.temperature_e.resize(N);
  results.temperature_i.resize(N);
  results.density.resize(N);
}

Loop::~Loop(void)
{
  //Destructor--free some stuff here
  delete heater;
  delete radiation_model;
}

void Loop::CalculateInitialConditions(void)
{
  int i = 0;
  int i_max = 100;
  double tol = 1e-2;
  double temperature_old = LARGEST_DOUBLE;
  double density_old = LARGEST_DOUBLE;
  double temperature,density;
  double radiative_loss;
  double error_temperature,error_density;
  double c1 = 2.0;
  double c2 = CalculateC2();
  double heat = heater->Get_Heating(0.0);

  while(i<i_max)
  {
    if(i > 0)
    {
      c1 = CalculateC1(temperature_old, temperature_old, density_old);
    }
    temperature = c2*pow(3.5*c1/(1.0 + c1)*pow(parameters.loop_length,2)*heat/(SPITZER_ELECTRON_CONDUCTIVITY + SPITZER_ION_CONDUCTIVITY),2.0/7.0);
    radiative_loss = radiation_model->GetPowerLawRad(log10(temperature));
    density = sqrt(heat/(radiative_loss*(1.0 + c1)));
    error_temperature = fabs(temperature - temperature_old)/temperature;
    error_density = fabs(density - density_old)/density;
    if(fmax(error_density,error_temperature) < tol)
    {
      break;
    }
    i++;
    temperature_old = temperature;
    density_old = density;
  }

  // Set current state in order pressure_e, pressure_i, density
  state.resize(3);
  state[0] = BOLTZMANN_CONSTANT*density*temperature;
  state[1] = parameters.boltzmann_correction*BOLTZMANN_CONSTANT*density*temperature;
  state[2] = density;

  //Save the results
  SaveResults(0,0.0);
}

void Loop::EvolveLoop(void)
{
  int i=1;
  double time = parameters.tau;
  double heat;

  while(time<parameters.total_time)
  {
    // Get current heating
    heat = heater->Get_Heating(time);
    // Solve Equations--update state
    // Save results
    SaveResults(i,time);
    //Update time and counter
    time += parameters.tau;
    i++;
  }
}

void Loop::PrintToFile(void)
{
  std::ofstream f;
  f.open(parameters.output_filename);
  for(int i=0;i<results.time.size();i++)
  {
    f << results.time[i] << "\t" << results.temperature_e[i] << "\t" << results.temperature_i[i] << "\t" << results.density[i] << "\t" << results.pressure_e[i] << "\t" << results.pressure_i[i] << "\t" << results.heat[i] << "\n";
  }
  f.close();
}

void Loop::SaveResults(int i,double time)
{
  // calculate parameters
  double temperature_e = state[0]/(BOLTZMANN_CONSTANT*state[2]);
  double temperature_i =  state[1]/(parameters.boltzmann_correction*BOLTZMANN_CONSTANT*state[2]);
  double heat = heater->Get_Heating(time);
  // Save results to results structure
  results.time[i] = time;
  results.heat[i] = heat;
  results.temperature_e[i] = temperature_e;
  results.temperature_i[i] = temperature_i;
  results.pressure_e[i] = state[0];
  results.pressure_i[i] = state[1];
  results.density[i] = state[2];
}

double Loop::CalculateC1(double temperature_e, double temperature_i, double density)
{
  double c1;
  double c1_eqm0 = 2.0;
  double c2 = CalculateC2();
  double grav_correction = 1.0;
  double loss_correction = 1.0;
  double scale_height = CalculateScaleHeight(temperature_e,temperature_i);
  double radiative_loss = radiation_model->GetPowerLawRad(log10(temperature_e));
  double density_eqm_2,density_ratio;

  if(parameters.use_c1_grav_correction)
  {
    grav_correction = exp(4.0*sin(_PI_/5.0)*parameters.loop_length/(_PI_*scale_height));
  }
  if(parameters.use_c1_loss_correction)
  {
    loss_correction = 1.95e-18*pow(temperature_e,-2.0/3.0)/radiative_loss;
  }

  density_eqm_2 = (SPITZER_ELECTRON_CONDUCTIVITY+SPITZER_ION_CONDUCTIVITY)*pow(temperature_e/c2,3.5)/(3.5*pow(parameters.loop_length,2)*c1_eqm0*loss_correction*grav_correction);
  density_ratio = pow(density,2)/density_eqm_2;

  if(density_ratio<1.0)
  {
    c1 = (2.0*c1_eqm0 + parameters.c1_cond0*(1.0/density_ratio - 1.0))/(1.0 + 1.0/density_ratio);
  }
  else
  {
    c1 = (2.0*c1_eqm0 + parameters.c1_rad0*(density_ratio - 1.0))/(1.0 + density_ratio);
  }

  return c1;
}

double Loop::CalculateC2(void)
{
  return 0.9;
}

double Loop::CalculateC3(void)
{
  return 0.6;
}

double Loop::CalculateScaleHeight(double temperature_e,double temperature_i)
{
  return BOLTZMANN_CONSTANT*(temperature_e + parameters.boltzmann_correction*temperature_i)/(parameters.ion_mass_correction*PROTON_MASS)/SOLAR_SURFACE_GRAVITY;
}

void Loop::CalculateAbundanceCorrection(double helium_to_hydrogen_ratio)
{
  double z_avg = (1.0 + 2.0*helium_to_hydrogen_ratio)/(1.0 + helium_to_hydrogen_ratio);
  parameters.boltzmann_correction = (1.0 + 1.0/z_avg)/2.0;
  parameters.ion_mass_correction = (1.0 + 4.0*helium_to_hydrogen_ratio)/(2.0 + 3.0*helium_to_hydrogen_ratio)*2.0*parameters.boltzmann_correction;
}
