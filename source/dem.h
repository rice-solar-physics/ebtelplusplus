/*
dem.h
Class definition for DEM object
*/

#ifndef DEM_H
#define DEM_H

#include "helper.h"
#include "loop.h"
#include "../Radiation_Model/source/radiation.h"
#include "../rsp_toolkit/source/xmlreader.h"
#include "../rsp_toolkit/source/file.h"
#include "../rsp_toolkit/source/constants.h"

class Dem{
private:
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

  // Calculate TR DEM
  //
  double CalculateDEMTR(int j,double density,double velocity,double pressure,double scale_height,double R_tr,double f_e);

public:
  // Default constructor
  //
  Dem(LOOP loop);

  // Destructor
  //
  ~Dem(void);

  // Calculate DEM
  //
  // Front end for DEM calculations; do any needed preprocessing
  // and calculations.
  //
  void CalculateDEM(int i);

  // Print results to file
  //
  void PrintToFile(int excess);
};

typedef Dem* DEM;

#endif
