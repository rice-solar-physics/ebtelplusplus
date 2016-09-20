/*
loop.h
Loop class definition
*/

#ifndef LOOP_H
#define LOOP_H

#include "helper.h"
#include "heater.h"
#include "../Radiation_Model/source/radiation.h"
#include "../rsp_toolkit/source/file.h"
#include "../rsp_toolkit/source/constants.h"

// Loop object
//
// Class for holding all of the information about the loop. It can
// be passed a configuration and after the initial conditions are
// set up, it is evolved through time.
//
class Loop {
private:

  /* Results structure */
  Results results;

  /* Instance of the <Heater> object */
  HEATER heater;

  /* Pointer to doc tree */
  tinyxml2::XMLDocument doc;

  /* Current state of the system */
  std::vector<double> __state;

  // Calculate c4
  //
  // Ratio of average to base velocity. Set to 1 for now
  //
  double CalculateC4(void);

  // Calculate coulomb collision frequency
  // @temperature_e electron temperature (in K)
  // @density number density (in cm^-3)
  //
  // Calculate the coulomb collision frequency for binary collisions
  // between electrons and ions according to Eq. 2.5e and Section 3 of
  // [Braginskii (1965)](http://adsabs.harvard.edu/abs/1965RvPP....1..205B).
  //
  // @return coulomb collision frequency (in s^-1)
  //
  double CalculateCollisionFrequency(double temperature_e,double density);

  // Calculate correction for He abundance
  //
  void CalculateAbundanceCorrection(double helium_to_hydrogen_ratio);

public:

  /* Instance of the <CRadiation> object */
  PRADIATION radiation_model;

  /* Parameter structure*/
  Parameters parameters;

  // Default constructor
  // @ebtel_config main configuration file
  // @rad_config configuration file for <radiation_model>; unused if using power-law radiative loss function
  //
  // Setup the loop object by reading in parameters from the configuration
  // file <ebtel_config> into the <parameters> structure. Additionally, the
  // <radiation_model> model object is created from the <rad_config> configuration
  // file. The constructor also creates the <heater> object for calculating
  // the heating profile.
  //
  Loop(char * ebtel_config,char * rad_config);

  /* Destructor */
  ~Loop(void);

  // Set initial conditions
  //
  void CalculateInitialConditions(void);

  // Print results to file
  //
  void PrintToFile(int excess);

  // Save results to structure
  //
  void SaveResults(int i, double time);

  // Return current state publicly
  //
  std::vector<double> GetState(void);

  // Set current state
  //
  void SetState(std::vector<double> state);

  // Calculate c1
  // @temperature_e electron temperature (in K)
  // @temperature_i ion temperature (in K)
  // @density number density (in cm^-3)
  //
  // Calculate the c1 parameter, the ratio between the
  // transition and coronal radiative losses
  //
  // @return c1 parameter
  //
  double CalculateC1(double temperature_e,double temperature_i,double density);

  // Calculate c2
  //
  double CalculateC2(void);

  // Calculate c3
  //
  double CalculateC3(void);

  // Calculate velocity
  //
  // Calculate the velocity using the base electron pressure and the enthalpy
  // flux as determined by our EBTEL equations.
  //
  double CalculateVelocity(double temperature_e,double temperature_i,double pressure_e);

  // Calculate temperature scale height
  //
  double CalculateScaleHeight(double temperature_e,double temperature_i);

  // Calculate thermal conduction
  //
  double CalculateThermalConduction(double temperature,double density,std::string species);

  // Calculate derivatives of EBTEL equations
  //
  std::vector<double> CalculateDerivs(std::vector<double> state,double time);
};

typedef Loop* LOOP;

#endif
