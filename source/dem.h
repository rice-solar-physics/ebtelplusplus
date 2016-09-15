/*
dem.h
Class definition for DEM object
*/

#ifndef DEM_H
#define DEM_H

#include "helper.h"
#include "../Radiation_Model/source/radiation.h"
#include "../rsp_toolkit/source/xmlreader.h"
#include "../rsp_toolkit/source/file.h"
#include "../rsp_toolkit/source/constants.h"

class Dem{
private:
  /* Ratio of average to apex temperature */
  double __c2;

  /* Ratio of base to apex temperature */
  double __c3;

  /* Loop half length */
  double __loop_length;

  /* Transition region DEM option */
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
  double CalculateDEMTR(void);

public:
  // Default constructor
  //
  Dem(tinyxml2::XMLElement * dem_node,PRADIATION radiation_model,size_t N,double loop_length,double c2,double c3);

  // Destructor
  //
  ~Dem(void);

  // Calculate DEM
  //
  // Front end for DEM calculations; do any needed preprocessing
  // and calculations.
  //
  void CalculateDEM(int i,double pressure,double density,double f_e,double c1);

  // Print results to file
  //
  void PrintToFile(std::string output_file,int excess);
};

typedef Dem* DEM;

#endif
