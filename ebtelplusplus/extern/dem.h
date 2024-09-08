/*
dem.h
Class definition for DEM object
*/

#ifndef DEM_H
#define DEM_H

#include "helper.h"
#include "loop.h"
#include "constants.h"

// DEM object
//
// Class for holding all of the methods needed to calculate
// the differential emission measure in the transition region
// and the corona. Requires the <loop> object for knowledge about
// the evolution of the coronal loop.
//
class Dem{
private:
  // Calculate TR DEM
  //
  double CalculateDEMTR(int j,double density,double velocity,double pressure,double scale_height,double R_tr,double f_e);

public:
  /* Loop object */
  LOOP loop;

  /* Method option for DEM TR calculation */
  bool use_new_method;

  /* Temperature range */
  std::vector<double> __temperature;

  /* Radiative loss */
  std::vector<double> __radiative_loss;

  /*Transition region DEM*/
  std::vector<std::vector<double> > dem_TR;

  /*Coronal DEM*/
  std::vector<std::vector<double> > dem_corona;

  // Default constructor
  //
  // Used when we don't want to actually do the DEM
  // calculation. Just a placeholder.
  //
  Dem(void);

  // Constructor
  // @loop <Loop> object that provides needed parameters and methods
  //
  // Setup Dem object to calculate differential emission measure in both the
  // transition region and the corona.
  //
  Dem(LOOP loop);

  // Destructor
  //
  ~Dem(void);

  // Calculate DEM
  // @i Timestep index
  //
  // Front end for DEM calculations. Calls methods to calculate both
  // the transition region and coronal DEM across the entire specified
  // temperature range.
  //
  void CalculateDEM(int i);

  // Print results to file
  // @num_steps number of steps taken by the integration routine
  //
  // Print coronal and transition region DEM arrays to separate files.
  // The filenames are the output filename as given in <loop>,
  // suffixed by `.dem_corona` and `.dem_tr`, respectively. The first
  // row of each file is the temperature vector, <__temperature>.
  //
  py::dict GetFinalResults(py::dict results, int num_steps);

};
// Pointer to the <Dem> class
typedef Dem* DEM;

#endif
