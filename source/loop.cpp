/*
loop.cpp
Loop object that will hold all information about the loop and be evolved in time.
*/

#include "loop.h"

Parameters Loop::parameters;
Terms Loop::terms;
HEATER Loop::heater;
PRADIATION Loop::radiation_model;

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
  parameters.adaptive_solver_error = std::stod(get_element_text(root,"adaptive_solver_error"));
  parameters.adaptive_solver_safety = std::stod(get_element_text(root,"adaptive_solver_safety"));
  parameters.saturation_limit = std::stod(get_element_text(root,"saturation_limit"));
  parameters.c1_cond0 = std::stod(get_element_text(root,"c1_cond0"));
  parameters.c1_rad0 = std::stod(get_element_text(root,"c1_rad0"));
  parameters.helium_to_hydrogen_ratio = std::stod(get_element_text(root,"helium_to_hydrogen_ratio"));
  //Boolean parameters
  parameters.force_single_fluid = string2bool(get_element_text(root,"force_single_fluid"));
  parameters.use_c1_loss_correction = string2bool(get_element_text(root,"use_c1_loss_correction"));
  parameters.use_c1_grav_correction = string2bool(get_element_text(root,"use_c1_grav_correction"));
  parameters.use_power_law_radiative_losses = string2bool(get_element_text(root,"use_power_law_radiative_losses"));
  parameters.use_flux_limiting = string2bool(get_element_text(root,"use_flux_limiting"));
  parameters.calculate_dem = string2bool(get_element_text(root,"calculate_dem"));
  parameters.use_adaptive_solver = string2bool(get_element_text(root,"use_adaptive_solver"));
  parameters.save_terms = string2bool(get_element_text(root,"save_terms"));
  //String parameters
  parameters.output_filename = get_element_text(root,"output_filename");

  //Estimate results array length
  parameters.N = int(std::ceil(parameters.total_time/parameters.tau))+1;

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

  //Initialize DEM object
  if(parameters.calculate_dem)
  {
    parameters.dem_options = get_element(root,"dem");
  }

  // Call the setup function
  Setup();
}

Loop::Loop(void)
{
  heater = new Heater();
}

Loop::~Loop(void)
{
  //Destructor--free some stuff here
  doc.Clear();
  delete heater;
  delete radiation_model;
}

void Loop::Setup(void)
{
  // Calculate needed He abundance corrections
  CalculateAbundanceCorrection(parameters.helium_to_hydrogen_ratio);

  //Reserve memory for results
  results.time.resize(parameters.N);
  results.heat.resize(parameters.N);
  results.pressure_e.resize(parameters.N);
  results.pressure_i.resize(parameters.N);
  results.temperature_e.resize(parameters.N);
  results.temperature_i.resize(parameters.N);
  results.density.resize(parameters.N);
  results.velocity.resize(parameters.N);
}

double Loop::GetMaxAllowedTimestep(void)
{
  return heater->minimum_duration*parameters.adaptive_solver_safety;
}

state_type Loop::GetState(void)
{
  return __state;
}

void Loop::SetState(state_type state)
{
  __state = state;
}

state_type Loop::CalculateInitialConditions(void)
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
  state_type state;

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
  state = {{ BOLTZMANN_CONSTANT*density*temperature, parameters.boltzmann_correction*BOLTZMANN_CONSTANT*density*temperature, density, temperature, temperature }};

  return state;
}

void Loop::PrintToFile(int num_steps)
{
  std::ofstream f;
  f.open(parameters.output_filename);
  for(int i=0;i<num_steps;i++)
  {
    f << results.time[i] << "\t" << results.temperature_e[i] << "\t" << results.temperature_i[i] << "\t" << results.density[i] << "\t" << results.pressure_e[i] << "\t" << results.pressure_i[i] << "\t" << results.velocity[i] << "\t" << results.heat[i] << "\n";
  }
  f.close();

  if(parameters.save_terms)
  {
    f.open(parameters.output_filename+".terms");
    for(int i=0;i<num_steps;i++)
    {
      f << terms.f_e[i] << "\t" << terms.f_i[i] << "\t" << terms.c1[i] << "\t" << terms.radiative_loss[i] << "\n";
    }
    f.close();
  }
}

void Loop::CalculateDerivs(const state_type &state, state_type &derivs, double time)
{
  double dpe_dt,dpi_dt,dn_dt,dTe_dt,dTi_dt;
  double psi_tr,psi_c,xi,R_tr,enthalpy_flux;

  double f_e = CalculateThermalConduction(state[3],state[2],"electron");
  double f_i = CalculateThermalConduction(state[4],state[2],"ion");
  double radiative_loss = radiation_model->GetPowerLawRad(std::log10(state[3]));
  double heat = heater->Get_Heating(time);
  double c1 = CalculateC1(state[3],state[4],state[2]);
  double c2 = CalculateC2();
  double c3 = CalculateC3();
  double collision_frequency = CalculateCollisionFrequency(state[3],state[2]);

  xi = state[0]/state[1];
  R_tr = c1*std::pow(state[2],2)*radiative_loss*parameters.loop_length;
  psi_tr = (f_e + R_tr - xi*f_i)/(1.0 + xi);
  psi_c = BOLTZMANN_CONSTANT*state[2]*collision_frequency*(state[4] - state[3]);
  enthalpy_flux = GAMMA_MINUS_ONE/GAMMA*(-f_e - R_tr + psi_tr);

  dpe_dt = GAMMA_MINUS_ONE*(heat*heater->partition + 1.0/parameters.loop_length*(psi_tr - R_tr*(1.0 + 1.0/c1))) + psi_c;
  dpi_dt = GAMMA_MINUS_ONE*(heat*(1.0 - heater->partition) - 1.0/parameters.loop_length*psi_tr) - psi_c;
  dn_dt = c2/(c3*parameters.loop_length*BOLTZMANN_CONSTANT*state[3])*enthalpy_flux;

  dTe_dt = state[3]*(1/state[0]*dpe_dt - 1/state[2]*dn_dt);
  dTi_dt = state[4]*(1/state[1]*dpi_dt - 1/state[2]*dn_dt);

  derivs[0] = dpe_dt;
  derivs[1] = dpi_dt;
  derivs[2] = dn_dt;
  derivs[3] = dTe_dt;
  derivs[4] = dTi_dt;
}

void Loop::SaveResults(int i,double time)
{
  // Get heating profile and velocity
  double heat = heater->Get_Heating(time);
  double velocity = CalculateVelocity(__state[3], __state[4], __state[0]);

  // Save results to results structure
  if(i >= parameters.N)
  {
    results.time.push_back(time);
    results.heat.push_back(heat);
    results.temperature_e.push_back(__state[3]);
    results.temperature_i.push_back(__state[4]);
    results.pressure_e.push_back(__state[0]);
    results.pressure_i.push_back(__state[1]);
    results.density.push_back(__state[2]);
    results.velocity.push_back(velocity);
  }
  else
  {
    results.time[i] = time;
    results.heat[i] = heat;
    results.temperature_e[i] = __state[3];
    results.temperature_i[i] = __state[4];
    results.pressure_e[i] = __state[0];
    results.pressure_i[i] = __state[1];
    results.density[i] = __state[2];
    results.velocity[i] = velocity;
  }
}

void Loop::SaveTerms(void)
{
  // Calculate terms
  double f_e = CalculateThermalConduction(__state[3], __state[2], "electron");
  double f_i = CalculateThermalConduction(__state[4], __state[2], "ion");
  double c1 = CalculateC1(__state[3], __state[4], __state[2]);
  double radiative_loss = radiation_model->GetPowerLawRad(std::log10(__state[3]));

  // Save terms
  terms.f_e.push_back(f_e);
  terms.f_i.push_back(f_i);
  terms.c1.push_back(c1);
  terms.radiative_loss.push_back(radiative_loss);
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

double Loop::CalculateCollisionFrequency(double temperature_e,double density)
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
