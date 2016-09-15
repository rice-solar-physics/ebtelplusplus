/*
dem.cpp
Methods for DEM class
*/

#include "dem.h"


Dem::Dem(tinyxml2::XMLElement * dem_node,PRADIATION radiation_model,size_t N,double loop_length,double c2,double c3)
{
  // Set some parameters for later calculations
  __loop_length = loop_length;
  __c2 = c2;
  __c3 = c3;
  use_new_method = string2bool(get_element_text(dem_node,"use_new_method"));
  // Configure temperature vector from inputs
  tinyxml2::XMLElement * temperature_node = get_element(dem_node,"temperature");
  int nbins = std::stoi(temperature_node->Attribute("bins"));
  double temperature_min = std::stod(temperature_node->Attribute("min"));
  double temperature_max = std::stod(temperature_node->Attribute("max"));
  double delta_temperature = (log10(temperature_max) - log10(temperature_min))/(nbins-1);
  // Set temperature bins and associated radiative losses
  __temperature.resize(nbins);
  __radiative_loss.resize(nbins);
  for(int i=0;i<nbins;i++)
  {
    __temperature[i] = temperature_min*pow(10.0,i*delta_temperature);
    __radiative_loss[i] = radiation_model->GetPowerLawRad(log10(__temperature[i]));
  }
  // Resize DEM arrays
  dem_TR.resize(N);
  dem_corona.resize(N);
  for(int i=0;i<N;i++)
  {
    dem_TR[i].resize(nbins);
    dem_corona[i].resize(nbins);
  }
}

Dem::~Dem(void)
{
  // Destructor--free some stuff here if needed
}

void Dem::CalculateDEM(int i,double temperature,double density,double f_e,double c1)
{
  double pressure = BOLTZMANN_CONSTANT*density*temperature;
  // Calculate coronal temperature range
  double temperature_corona_max = fmax(temperature/__c2,1.1e+4);
  double temperature_corona_min = fmax(temperature*(2.0 - 1.0/__c2),1.0e+4);
  // Calculate coronal emission
  double delta_temperature = pow(10.0,0.5)*temperature_corona_max - pow(10.0,-0.5)*temperature_corona_min;
  double coronal_emission = 2.0*pow(density,2)*__loop_length/delta_temperature;

  for(int j=0;j<__temperature.size();j++)
  {
    // Coronal DEM
    dem_corona[i][j] = 0.0;
    if(__temperature[j]<=temperature_corona_max && __temperature[j]>=temperature_corona_min)
    {
      dem_corona[i][j] = coronal_emission;
    }
    // Transition Region
    dem_TR[i][j] = 0.0;
    if(__temperature[j]<__c3/__c2*temperature)
    {
      dem_TR[i][j] = CalculateDEMTR();
    }
  }
}

void Dem::PrintToFile(std::string output_file,int excess)
{
  // Trim zeros
  for(int i=0;i<excess;i++)
  {
    dem_TR.pop_back();
    dem_corona.pop_back();
  }

  // Open file streams
  std::ofstream f_corona;
  std::ofstream f_tr;
  f_corona.open(output_file+".dem_corona");
  f_tr.open(output_file+".dem_tr");
  // First row of each file is the temperature array
  for(int j=0;j<__temperature.size();j++)
  {
    f_corona << __temperature[j] << "\t";
    f_tr << __temperature[j] << "\t";
  }
  f_corona << "\n";
  f_tr << "\n";
  // Print TR and corona DEM at each timestep
  for(int i=0;i<dem_TR.size();i++)
  {
    for(int j=0;j<dem_TR[i].size();j++)
    {
      f_corona << dem_corona[i][j] << "\t";
      f_tr << dem_TR[i][j] << "\t";
    }
    f_corona << "\n";
    f_tr << "\n";
  }

  // Close the file
  f_corona.close();
  f_tr.close();
}

double Dem::CalculateDEMTR(void)
{
  return 1.0;
}
