/*
ebtel++
A code for computing the evolution of dynamically heated, spatially-averaged solar coronal loops.
*/

#include <malloc.h>
#include <time.h>
#include "boost/program_options.hpp"
#include "loop.h"

int main(int argc, char *argv[])
{

  //Declarations
  char rad_config[256];
  char ebtel_config[256];
  LOOP loop;

  //Parse command line options with boost
  namespace po = boost::program_options;
  po::options_description description("ebtel++ A code for efficiently computing the evolution of dynamically-heated coronal loops. Based on the model of Klimchuk et al. (2008).");
  description.add_options()
    ("help,h","The help message")
    ("quiet,q",po::bool_switch()->default_value(false),"Suppress output.")
    ("ebtel_config,c",po::value<std::string>()->default_value("config/ebtel.example.cfg.xml"),"Configuration file for EBTEL.")
    ("rad_config,r",po::value<std::string>()->default_value("config/radiation.example.cfg.xml"),"Configuration file for radiation class");
  po::variables_map vm;
  po::store(po::command_line_parser(argc,argv).options(description).run(), vm);
  if(vm.count("help"))
  {
  	std::cout << description;
  	return 0;
  }
  po::notify(vm);

  //Copy parameter to char array
  std::strcpy(rad_config,vm["rad_config"].as<std::string>().c_str());
  std::strcpy(ebtel_config,vm["ebtel_config"].as<std::string>().c_str());

  //Create loop object
  loop = new Loop(ebtel_config,rad_config);

  //Set initional conditions of the loop
  loop->CalculateInitialConditions();
  //Evolve loop
  loop->EvolveLoop();
  //Print results to file
  loop->PrintToFile();

  //Destroy loop object
  delete loop;

  return 0;
}
