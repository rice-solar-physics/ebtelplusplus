/*
helper.h
General purpose includes to be used everywhere
*/

#ifndef HELPER_H
#define HELPER_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>

struct Parameters {
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
};

struct Results {
  std::vector<double> time;
  std::vector<double> heat;
};

#endif
