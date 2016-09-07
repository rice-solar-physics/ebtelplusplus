/*
loop.cpp
Loop object that will hold all information about the loop and be evolved in time.
*/

#include "loop.h"

Loop::Loop(char *ebtel_config, PRADIATION radiation_model)
{
  tinyxml2::XMLDocument doc;

  //Open file
  tinyxml2::XMLError load_ok = doc.LoadFile(ebtel_config);
  if(load_ok != 0)
  {
  	printf("Failed to load XML configuration file %s.\n",ebtel_config);
  	//TODO: Exit or break out from here
  }
  //Parse file and read into data structure
  tinyxml2::XMLElement *root = doc.FirstChildElement();
  //Numeric parameters
  parameters.total_time = std::stod(get_element_text(root,"total_time"));
  parameters.tau = std::stod(get_element_text(root,"tau"));
  parameters.loop_length = std::stod(get_element_text(root,"loop_length"));
  parameters.rka_error = std::stod(get_element_text(root,"loop_length"));
  //Boolean parameters
  parameters.use_c1_loss_correction = string2bool(get_element_text(root,"use_c1_loss_correction"));
  parameters.use_c1_grav_correction = string2bool(get_element_text(root,"use_c1_grav_correction"));
  parameters.use_power_law_radiative_losses = string2bool(get_element_text(root,"use_power_law_radiative_losses"));
  parameters.use_spitzer_conductivity = string2bool(get_element_text(root,"use_spitzer_conductivity"));
  //String parameters
  parameters.solver = get_element_text(root,"solver");

  //Initialize heating object
  //Initialize DEM object
  //TODO: implement after everything else is in place

  doc.Clear();
}

Loop::~Loop(void)
{
  //Destructor--free some stuff here
}

void Loop::CalculateInitialConditions(void)
{
  //Calculate initial conditions
  //Set initial state
}
