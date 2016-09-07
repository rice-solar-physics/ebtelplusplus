/*
heater.h
Class definition for the heating object
*/

#ifndef HEATER_H
#define HEATER_H

#include "helper.h"
#include "../rsp_toolkit/source/xmlreader.h"

class Heater {
private:

  /*Background heating rate*/
  double background;

  /* Number of events */
  int num_events;

  /*Starting time of the rise phase*/
  std::vector<double> time_start_rise;

  /*Ending time of the rise phase*/
  std::vector<double> time_end_rise;

  /*Starting time of the decay phase*/
  std::vector<double> time_start_decay;

  /*Ending time of the decay phase*/
  std::vector<double> time_end_decay;

  /*Magnitudes of the events*/
  std::vector<double> magnitude;

public:

  /* Partition of energy between electrons and ions; 1 is electron heating, 0 is ion heating */
  double partition;

  // Default constructor
  // @heating_node <tinyxml2::XMLElement> object holding the heating information
  //
  Heater(tinyxml2::XMLElement * heating_node);

  /* Destructor */
  ~Heater(void);

  // Get heating at time <t>
  // @t current time
  //
  // Given the heating profile specified in the configuration file,
  // return the heating rate at the given time <t>
  //
  // @return heating rate at time t.
  //
  double Get_Heating(double t);

};
typedef Heater* HEATER;

#endif
