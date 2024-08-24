/*
heater.cpp
Class defnition for heater object
*/

#include "heater.h"

Heater::Heater(py::dict heating_config)
{
  background = heating_config["background"].cast<float>();
  partition = heating_config["partition"].cast<float>();
  py::list event_list = heating_config["events"];
  for(auto event:event_list)
  {
    time_start_rise.push_back(event["event"]["rise_start"].cast<float>());
    time_end_rise.push_back(event["event"]["rise_end"].cast<float>());
    time_start_decay.push_back(event["event"]["decay_start"].cast<float>());
    time_end_decay.push_back(event["event"]["decay_end"].cast<float>());
    magnitude.push_back(event["event"]["magnitude"].cast<float>());
  }
  num_events = magnitude.size();
}

Heater::Heater(void)
{
  // Default constructor
}

Heater::~Heater(void)
{
  //Destructor--free some stuff here
}

double Heater::Get_Heating(double time)
{
  double heat = background;
  for(int i=0;i<num_events;i++)
  {
    if(time >= time_start_rise[i] && time < time_end_rise[i])
    {
      heat += magnitude[i]*(time - time_start_rise[i])/(time_end_rise[i] - time_start_rise[i]);
    }
    else if(time >= time_end_rise[i] && time < time_start_decay[i])
    {
      heat += magnitude[i];
    }
    else if(time >= time_start_decay[i] && time < time_end_decay[i])
    {
      heat += magnitude[i]*(time_end_decay[i] - time)/(time_end_decay[i] - time_start_decay[i]);
    }
  }

  return heat;
}

double Heater::Get_Time_To_Next_Heating_Change(double time)
{
  double tau = std::numeric_limits<double>::max();
  for(int i=0;i<num_events;i++)
  {
    if(time < time_start_rise[i])
    {
      tau = std::fmin( tau, time_start_rise[i] - time);
    }
    else if(time < time_end_rise[i])
    {
      tau = std::fmin( tau, time_end_rise[i] - time);
    }
    else if(time < time_start_decay[i])
    {
      tau = std::fmin( tau, time_start_decay[i] - time);
    }
    else if(time < time_end_decay[i])
    {
      tau = std::fmin( tau, time_end_decay[i] - time);
    }
  }

  return tau;
}
