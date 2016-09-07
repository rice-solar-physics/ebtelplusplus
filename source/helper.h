/*
helper.h
General purpose includes to be used everywhere
*/

#ifndef HELPER_H
#define HELPER_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>

struct Parameters {
  double total_time,tau;
  double loop_length;
  double rka_error;
  bool use_c1_loss_correction,use_c1_grav_correction;
  bool use_power_law_radiative_losses,use_spitzer_conductivity;
  std::string solver;
};

#endif
