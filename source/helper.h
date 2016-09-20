/*
helper.h
General purpose includes to be used everywhere
*/

#ifndef HELPER_H
#define HELPER_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../rsp_toolkit/source/xmlreader.h"

struct Parameters {
  /* Total simulation time (in s) */
  double total_time;
  /* Timestep (in s); when using adaptive solver, the initial timestep*/
  double tau;
  /* Loop half length (in cm) */
  double loop_length;
  /* Truncation error tolerance for adaptive solver */
  double rka_error;
  /* Heat flux saturation limit; 1/6 is a typical value */
  double saturation_limit;
  /* Nominal conductive c1 value; 6.0 is recommended*/
  double c1_cond0;
  /* Nominal radiative c1 value; 0.6 is recommended */
  double c1_rad0;
  /* Switch for radiative loss correction to c1 */
  bool use_c1_loss_correction;
  /* Switch for gravitational correction to c1 */
  bool use_c1_grav_correction;
  /* Switch for power-law radiative loss function */
  bool use_power_law_radiative_losses;
  /* Switch for classical Spitzer conductivity */
  bool use_spitzer_conductivity;
  /* Switch for calculating DEM; if True, runtimes will be much longer */
  bool calculate_dem;
  /* Solver option; use `euler`, `rk4`, or `rka4` */
  std::string solver;
  /* Path to output file */
  std::string output_filename;
  /* XML node holding DEM calculation parameters */
  tinyxml2::XMLElement * dem_options;
  /* Correction to ion equation of state */
  double boltzmann_correction;
  /* Ion mass correction to account for He abundance */
  double ion_mass_correction;
  /* Number of grid points */
  size_t N;
};

struct Results {
  /* Time (in s) */
  std::vector<double> time;
  /* Electron temperature (in K) */
  std::vector<double> temperature_e;
  /* Ion temperature (in K) */
  std::vector<double> temperature_i;
  /* Electron pressure (in dyne cm^-2 s^-1) */
  std::vector<double> pressure_e;
  /* Ion pressure (in dyne cm^-2 s^-1) */
  std::vector<double> pressure_i;
  /* Number density (in cm^-3) */
  std::vector<double> density;
  /* Heating rate (in erg cm^-3 s^-1) */
  std::vector<double> heat;
};

#endif
