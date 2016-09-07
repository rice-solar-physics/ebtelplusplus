/*
loop.h
Loop class definition
*/

#ifndef LOOP_H
#define LOOP_H

#include "helper.h"
#include "../Radiation_Model/source/radiation.h"
#include "../rsp_toolkit/source/file.h"

// Loop class
//
// Class for holding all of the information about the loop. It can
// be passed a configuration and after the initial conditions are
// set up, it is evolved through time.
//
class Loop {
private:

  // Parameter structure
  Parameters parameters;

public:

  // Default constructor
  // @ebtel_config name of main configuration file
  // @radiation_model instance of <CRadiation> class
  //
  Loop(char * ebtel_config,PRADIATION radiation_model);

  /* Destructor */
  ~Loop(void);

  // Set initial conditions
  //
  void CalculateInitialConditions(void);
};

typedef Loop* LOOP;

#endif
