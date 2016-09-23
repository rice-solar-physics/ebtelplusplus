/*
loop.cpp
Loop object that will hold all information about the loop and be evolved in time.
*/

#include "loop.h"

Loop::Loop(char *ebtel_config, char *rad_config)
{
  tinyxml2::XMLElement *root;
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
  parameters.rka_error = std::stod(get_element_text(root,"rka_error"));
  parameters.saturation_limit = std::stod(get_element_text(root,"saturation_limit"));
  parameters.c1_cond0 = std::stod(get_element_text(root,"c1_cond0"));
  parameters.c1_rad0 = std::stod(get_element_text(root,"c1_rad0"));
  helium_to_hydrogen_ratio = std::stod(get_element_text(root,"helium_to_hydrogen_ratio"));
  //Boolean parameters
  parameters.force_single_fluid = string2bool(get_element_text(root,"force_single_fluid"));
  parameters.use_c1_loss_correction = string2bool(get_element_text(root,"use_c1_loss_correction"));
  parameters.use_c1_grav_correction = string2bool(get_element_text(root,"use_c1_grav_correction"));
  parameters.use_power_law_radiative_losses = string2bool(get_element_text(root,"use_power_law_radiative_losses"));
  parameters.use_flux_limiting = string2bool(get_element_text(root,"use_flux_limiting"));
  parameters.calculate_dem = string2bool(get_element_text(root,"calculate_dem"));
  //String parameters
  parameters.solver = get_element_text(root,"solver");
  parameters.output_filename = get_element_text(root,"output_filename");

  //Estimate results array length
  parameters.N = int(std::ceil(parameters.total_time/parameters.tau));

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
  heater = new Heater(get_element(root,"heating"));

  // Calculate needed He abundance corrections
  CalculateAbundanceCorrection(helium_to_hydrogen_ratio);

  //Initialize DEM object
  if(parameters.calculate_dem)
  {
    parameters.dem_options = get_element(root,"dem");
  }

  //Reserve memory for results
  results.time.resize(parameters.N);
  results.heat.resize(parameters.N);
  results.pressure_e.resize(parameters.N);
  results.pressure_i.resize(parameters.N);
  results.temperature_e.resize(parameters.N);
  results.temperature_i.resize(parameters.N);
  results.density.resize(parameters.N);
}

Loop::~Loop(void)
{
  //Destructor--free some stuff here
  doc.Clear();
  delete heater;
  delete radiation_model;
}

std::vector<double> Loop::GetState(void)
{
  return __state;
}

void Loop::SetState(std::vector<double> state)
{
  __state = state;
}

void Loop::CalculateInitialConditions(void)
{
  int i = 0;
  int i_max = 100;
  double tol = 1e-2;
  double temperature_old = (double)LARGEST_DOUBLE;
  double density_old = (double)LARGEST_DOUBLE;
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
    temperature = c2*std::pow(3.5*c1/(1.0 + c1)*std::pow(parameters.loop_length,2)*heat/(SPITZER_ELECTRON_CONDUCTIVITY + SPITZER_ION_CONDUCTIVITY),2.0/7.0);
    radiative_loss = radiation_model->GetPowerLawRad(std::log10(temperature));
    density = std::sqrt(heat/(radiative_loss*(1.0 + c1)));
    error_temperature = std::abs(temperature - temperature_old)/temperature;
    error_density = std::abs(density - density_old)/density;
    if(std::fmax(error_density,error_temperature) < tol)
    {
      break;
    }
    i++;
    temperature_old = temperature;
    density_old = density;
  }

  // Set current state in order pressure_e, pressure_i, density
  __state.resize(3);
  __state[0] = BOLTZMANN_CONSTANT*density*temperature;
  __state[1] = parameters.boltzmann_correction*BOLTZMANN_CONSTANT*density*temperature;
  __state[2] = density;

  //Save the results
  SaveResults(0,0.0);
}

void Loop::PrintToFile(int excess)
{
  int i;
  // Trim zeroes
  for(i=0;i<excess;i++)
  {
    results.time.pop_back();
    results.temperature_e.pop_back();
    results.temperature_i.pop_back();
    results.density.pop_back();
    results.pressure_e.pop_back();
    results.pressure_i.pop_back();
    results.heat.pop_back();
  }

  std::ofstream f;
  f.open(parameters.output_filename);
  for(i=0;i<results.time.size();i++)
  {
    f << results.time[i] << "\t" << results.temperature_e[i] << "\t" << results.temperature_i[i] << "\t" << results.density[i] << "\t" << results.pressure_e[i] << "\t" << results.pressure_i[i] << "\t" << results.heat[i] << "\n";
  }
  f.close();
}

std::vector<double> Loop::CalculateDerivs(std::vector<double> state,double time)
{
  std::vector<double> derivs(3);
  double dpe_dt,dpi_dt,dn_dt;
  long double psi_tr,psi_c,xi,R_tr,enthalpy_flux;

  long double temperature_e = state[0]/(BOLTZMANN_CONSTANT*state[2]);
  long double temperature_i = state[1]/(BOLTZMANN_CONSTANT*parameters.boltzmann_correction*state[2]);

  double f_e = CalculateThermalConduction(temperature_e,state[2],"electron");
  double f_i = CalculateThermalConduction(temperature_i,state[2],"ion");
  double radiative_loss = radiation_model->GetPowerLawRad(std::log10(temperature_e));
  double heat = heater->Get_Heating(time);
  double c1 = CalculateC1(temperature_e,temperature_i,state[2]);
  double c2 = CalculateC2();
  double c3 = CalculateC3();
  double collision_frequency = CalculateCollisionFrequency(temperature_e,state[2]);

  xi = state[0]/state[1];
  R_tr = c1*std::pow(state[2],2)*radiative_loss*parameters.loop_length;
  psi_tr = (f_e + R_tr - xi*f_i)/(1.0 + xi);
  psi_c = BOLTZMANN_CONSTANT*state[2]*collision_frequency*(temperature_i - temperature_e);
  enthalpy_flux = GAMMA_MINUS_ONE/GAMMA*(-f_e - R_tr + psi_tr);

  dpe_dt = GAMMA_MINUS_ONE*(heat*heater->partition + 1.0/parameters.loop_length*(psi_tr - R_tr*(1.0 + 1.0/c1))) + psi_c;
  dpi_dt = GAMMA_MINUS_ONE*(heat*(1.0 - heater->partition) - 1.0/parameters.loop_length*psi_tr) - psi_c;
  dn_dt = c2/(c3*parameters.loop_length*BOLTZMANN_CONSTANT*temperature_e)*enthalpy_flux;

  derivs[0] = dpe_dt;
  derivs[1] = dpi_dt;
  derivs[2] = dn_dt;

  return derivs;
}

void Loop::SaveResults(int i,double time)
{
  // calculate parameters
  double heat = heater->Get_Heating(time);
  double temperature_e = __state[0]/(BOLTZMANN_CONSTANT*__state[2]);
  double temperature_i = __state[1]/(BOLTZMANN_CONSTANT*parameters.boltzmann_correction*__state[2]);

  // Save results to results structure
  if(i >= parameters.N)
  {
    results.time.push_back(time);
    results.heat.push_back(heat);
    results.temperature_e.push_back(temperature_e);
    results.temperature_i.push_back(temperature_i);
    results.pressure_e.push_back(__state[0]);
    results.pressure_i.push_back(__state[1]);
    results.density.push_back(__state[2]);
  }
  else
  {
    results.time[i] = time;
    results.heat[i] = heat;
    results.temperature_e[i] = temperature_e;
    results.temperature_i[i] = temperature_i;
    results.pressure_e[i] = __state[0];
    results.pressure_i[i] = __state[1];
    results.density[i] = __state[2];
  }
}

double Loop::CalculateThermalConduction(double temperature, double density, std::string species)
{
  double kappa,mass,k_B;
  double f_c,f;
  double c2 = CalculateC2();

  if(species.compare("electron")==0)
  {
    kappa = SPITZER_ELECTRON_CONDUCTIVITY;
    mass = ELECTRON_MASS;
    k_B = BOLTZMANN_CONSTANT;
  }
  else
  {
    kappa = SPITZER_ION_CONDUCTIVITY;
    mass = parameters.ion_mass_correction*PROTON_MASS;
    k_B = parameters.boltzmann_correction*BOLTZMANN_CONSTANT;
  }

  f_c = -2.0/7.0*kappa*std::pow(temperature/c2,3.5)/parameters.loop_length;

  if(parameters.use_flux_limiting)
  {
    double f_s = -parameters.saturation_limit*1.5/std::sqrt(mass)*density*std::pow(k_B*temperature,1.5);
    f = -f_c*f_s/std::sqrt(std::pow(f_c,2) + std::pow(f_s,2));
  }
  else
  {
    f = f_c;
  }

  return f;
}

double Loop::CalculateCollisionFrequency(double temperature_e, double temperature_i,double density)
{
  if(parameters.force_single_fluid)
  {
    // TODO: explain why this value works
    return 0.9;
  }
  else
  {
    // TODO: find a reference for this formula
    double coulomb_logarithm = 23.0 - std::log(std::sqrt(density/1.0e+13)*std::pow(BOLTZMANN_CONSTANT*temperature_e/(1.602e-9),-1.5));
    return 16.0*SQRT_PI/3.0*ELECTRON_CHARGE_POWER_4/(parameters.ion_mass_correction*PROTON_MASS*ELECTRON_MASS)*std::pow(2.0*BOLTZMANN_CONSTANT*temperature_e/ELECTRON_MASS,-1.5)*density*coulomb_logarithm;
  }
}

double Loop::CalculateC1(double temperature_e, double temperature_i, double density)
{
  double c1;
  double density_eqm_2,density_ratio;

  double c1_eqm0 = 2.0;
  double c2 = CalculateC2();
  double grav_correction = 1.0;
  double loss_correction = 1.0;
  double scale_height = CalculateScaleHeight(temperature_e,temperature_i);
  double radiative_loss = radiation_model->GetPowerLawRad(std::log10(temperature_e));

  if(parameters.use_c1_grav_correction)
  {
    grav_correction = std::exp(4.0*std::sin(_PI_/5.0)*parameters.loop_length/(_PI_*scale_height));
  }
  if(parameters.use_c1_loss_correction)
  {
    loss_correction = 1.95e-18/std::pow(temperature_e,2.0/3.0)/radiative_loss;
  }

  density_eqm_2 = (SPITZER_ELECTRON_CONDUCTIVITY + SPITZER_ION_CONDUCTIVITY)*std::pow(temperature_e/c2,3.5)/(3.5*std::pow(parameters.loop_length,2)*c1_eqm0*loss_correction*grav_correction*radiative_loss);
  density_ratio = std::pow(density,2)/density_eqm_2;

  if(density_ratio<1.0)
  {
    c1 = (2.0*c1_eqm0 + parameters.c1_cond0*(1.0/density_ratio - 1.0))/(1.0 + 1.0/density_ratio);
  }
  else
  {
    c1 = (2.0*c1_eqm0 + parameters.c1_rad0*(density_ratio - 1.0))/(1.0 + density_ratio);
  }

  return c1*loss_correction*grav_correction;
}

double Loop::CalculateC2(void)
{
  return 0.9;
}

double Loop::CalculateC3(void)
{
  return 0.6;
}

double Loop::CalculateC4(void)
{
  return 1.0;
}

double Loop::CalculateScaleHeight(double temperature_e,double temperature_i)
{
  return BOLTZMANN_CONSTANT*(temperature_e + parameters.boltzmann_correction*temperature_i)/(parameters.ion_mass_correction*PROTON_MASS)/((double)SOLAR_SURFACE_GRAVITY);
}

void Loop::CalculateAbundanceCorrection(double helium_to_hydrogen_ratio)
{
  double z_avg = (1.0 + 2.0*helium_to_hydrogen_ratio)/(1.0 + helium_to_hydrogen_ratio);
  parameters.boltzmann_correction = (1.0 + 1.0/z_avg)/2.0;
  parameters.ion_mass_correction = (1.0 + 4.0*helium_to_hydrogen_ratio)/(2.0 + 3.0*helium_to_hydrogen_ratio)*2.0*parameters.boltzmann_correction;
}

double Loop::CalculateVelocity(double temperature_e, double temperature_i, double pressure_e)
{
  double c4 = CalculateC4();
  double density = pressure_e/(BOLTZMANN_CONSTANT*temperature_e);
  double c1 = CalculateC1(temperature_e,temperature_i,density);
  double R_tr = c1*std::pow(density,2)*radiation_model->GetPowerLawRad(std::log10(temperature_e))*parameters.loop_length;
  double fe = CalculateThermalConduction(temperature_e,density,"electron");
  double fi = CalculateThermalConduction(temperature_i,density,"ion");
  double sc = CalculateScaleHeight(temperature_e,temperature_i);
  double xi = temperature_e/temperature_i/parameters.boltzmann_correction;

  double coefficient = c4*xi*GAMMA_MINUS_ONE/(GAMMA*(xi+1));
  double pressure_e_0 = pressure_e*std::exp(2.0*parameters.loop_length*std::sin(_PI_/5.0)/(_PI_*sc));
  double enthalpy_flux = -(fe + fi + R_tr);

  return coefficient*enthalpy_flux/pressure_e_0;
}
