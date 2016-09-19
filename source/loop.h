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

// Loop class
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
  //
  double CalculateCollisionFrequency(double temperature_e,double density);

  // Calculate correction for He abundance
  //
  void CalculateAbundanceCorrection(double helium_to_hydrogen_ratio);

  // Calculate derivatives of EBTEL equations
  //
  std::vector<double> CalculateDerivs(std::vector<double> state,double time);

public:

  /* Instance of the <CRadiation> object */
  PRADIATION radiation_model;

  /* Parameter structure*/
  Parameters parameters;

  // Default constructor
  // @ebtel_config name of main configuration file
  // @radiation_model instance of <CRadiation> class
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

  // Euler solver
  //
  std::vector<double> EulerSolver(std::vector<double> state, double time, double tau);

  // Fourth-order Runge-Kutta solver
  //
  std::vector<double> RK4Solver(std::vector<double> state, double time, double tau);

  // Adaptive time-stepper for RK4 solver
  //
  std::vector<double> RKA4Solver(std::vector<double> state, double time, double tau);
};

typedef Loop* LOOP;

#endif
