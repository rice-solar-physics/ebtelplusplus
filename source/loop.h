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

  /* Parameter structure*/
  Parameters parameters;

  /* Results structure */
  Results results;

  /* Instance of the <Heater> object */
  HEATER heater;

  /* Instance of the <CRadiation> object */
  PRADIATION radiation_model;

  /* Current state of the system */
  std::vector<double> state;

  /* Number of entries in results */
  size_t N;

  // Save results to structure
  //
  void SaveResults(int i, double time);

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

  // Calculate temperature scale height
  //
  double CalculateScaleHeight(double temperature_e,double temperature_i);

  // Calculate correction for He abundance
  //
  void CalculateAbundanceCorrection(double helium_to_hydrogen_ratio);

public:

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

  // Evolve loop in time
  //
  void EvolveLoop(void);

  // Print results to file
  //
  void PrintToFile(void);
};

typedef Loop* LOOP;

#endif
