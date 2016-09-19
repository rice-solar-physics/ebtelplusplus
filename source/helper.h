/*
helper.h
General purpose includes to be used everywhere
*/

#ifndef HELPER_H
#define HELPER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../rsp_toolkit/source/xmlreader.h"

struct Parameters {
  // Input parameters
  double total_time,tau;
  double loop_length;
  double rka_error;
  double saturation_limit;
  double c1_cond0,c1_rad0;
  bool use_c1_loss_correction,use_c1_grav_correction;
  bool use_power_law_radiative_losses,use_spitzer_conductivity;
  bool calculate_dem;
  std::string solver;
  std::string output_filename;
  tinyxml2::XMLElement * dem_options;
  // Calculated parameters
  double boltzmann_correction;
  double ion_mass_correction;
  size_t N;
};

struct Results {
  std::vector<double> time;
  std::vector<double> temperature_e;
  std::vector<double> temperature_i;
  std::vector<double> pressure_e;
  std::vector<double> pressure_i;
  std::vector<double> density;
  std::vector<double> heat;
};

#endif
