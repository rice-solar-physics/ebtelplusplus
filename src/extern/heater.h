/*
heater.h
Class definition for the heating object
*/

#ifndef HEATER_H
#define HEATER_H

#include "helper.h"
#include "constants.h"

// Heater object
//
// Class for configuring time-dependent heating profiles.
// Accepts a properly formatted XML node and
// calculates the heating rate at any time. Heating profiles
// must be specified in terms of <num_events> heating pulses
// plus a static background <background>. You can also initialize a blank
// object and set the event parameters later on.
//
class Heater {

public:

  /*Background heating rate (in erg cm^-3 s^-1) */
  double background;

  /* Number of events */
  int num_events;

  /*Starting time of the rise phase (in s) */
  std::vector<double> time_start_rise;

  /*Ending time of the rise phase (in s) */
  std::vector<double> time_end_rise;

  /*Starting time of the decay phase (in s) */
  std::vector<double> time_start_decay;

  /*Ending time of the decay phase (in s) */
  std::vector<double> time_end_decay;

  /*Heating rates of the events (in erg cm^-3 s^-1) */
  std::vector<double> rate;

  /* Partition of energy between electrons and ions; 1 corresponds to pure electron heating and 0 pure ion heating. For a single-fluid treatment, use 0.5 */
  double partition;

  // Constructor
  // @heating_node XML node holding the heating information
  //
  Heater(py::dict heating_config);

  // Default constructor
  //
  Heater(void);

  /* Destructor */
  ~Heater(void);

  // Get heating at time <time>
  // @time current time (in s)
  //
  // Given the heating profile specified by the configuration file,
  // return the heating rate at the given time <time>
  //
  // @return heating rate at time <time> (in erg cm^-3 s^-1)
  //
  double Get_Heating(double time);

  // Get the time until the next change in heating profile <time>
  // @time current time (in s)
  //
  // Calculates the time until the next heating_start_rise, heating_end_rise,
  // heating_start_decay, or heating_end_decay. Used to ensure the time stepper
  // doesn't skip over any transitions in the heating function.
  //
  // @return time until next change in the loop heating function (in s)
  //
  double Get_Time_To_Next_Heating_Change(double time);
};
// Pointer to the <Heater> class
typedef Heater* HEATER;

#endif
