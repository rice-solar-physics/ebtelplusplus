/*
heater.cpp
Class defnition for heater object
*/

#include "heater.h"

Heater::Heater(tinyxml2::XMLElement * heating_node)
{
  //Set basic parameters
  background = std::stod(get_element_text(heating_node,"background"));
  partition = std::stod(get_element_text(heating_node,"partition"));

  //Set heating parameters
  tinyxml2::XMLElement * events = get_element(heating_node,"events");
  for(tinyxml2::XMLElement *child = events->FirstChildElement();child != NULL;child=child->NextSiblingElement())
  {
    time_start_rise.push_back(std::stod(child->Attribute("rise_start")));
    time_end_rise.push_back(std::stod(child->Attribute("rise_end")));
    time_start_decay.push_back(std::stod(child->Attribute("decay_start")));
    time_end_decay.push_back(std::stod(child->Attribute("decay_end")));
    magnitude.push_back(std::stod(child->Attribute("magnitude")));
  }
  num_events = magnitude.size();

}

Heater::~Heater(void)
{
  //Destructor--free some stuff here
}

double Heater::Get_Heating(double t)
{
  double heat = background;
  for(int i=0;i<num_events;i++)
  {
    if(t >= time_start_rise[i] && t < time_end_rise[i])
    {
      heat += magnitude[i]*(t - time_start_rise[i])/(time_end_rise[i] - time_start_rise[i]);
    }
    else if(t >= time_end_rise[i] && t < time_start_decay[i])
    {
      heat += magnitude[i];
    }
    else if(t >= time_start_decay[i] && t < time_end_decay[i])
    {
      heat += magnitude[i]*(time_end_decay[i] - t)/(time_end_decay[i] - time_start_decay[i]);
    }
  }

  return heat;
}