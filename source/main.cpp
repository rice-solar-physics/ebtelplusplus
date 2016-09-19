/*
ebtel++
A code for computing the evolution of dynamically heated, spatially-averaged solar coronal loops.
*/

#include <malloc.h>
#include <time.h>
#include "boost/program_options.hpp"
#include "loop.h"
#include "dem.h"

int main(int argc, char *argv[])
{

  //Declarations
  int i;
  double time,tau;
  std::vector<double> state;
  char rad_config[256],ebtel_config[256];
  LOOP loop;
  DEM dem;

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

  // Create DEM object if needed
  if(loop->parameters.calculate_dem)
  {
    dem = new Dem(loop);
  }

  //Evolve loop
  time = loop->parameters.tau;
  tau = loop->parameters.tau;
  i = 1;
  state = loop->GetState();
  while(time<loop->parameters.total_time)
  {
    // Solve Equations--update state
    if(loop->parameters.solver.compare("euler")==0)
    {
      state = loop->EulerSolver(state,time,tau);
    }
    else if(loop->parameters.solver.compare("rk4")==0)
    {
      state = loop->RK4Solver(state,time,tau);
    }
    else if(loop->parameters.solver.compare("rka4")==0)
    {
      std::vector<double> _tmp_state;
      _tmp_state = loop->RKA4Solver(state,time,tau);
      tau = _tmp_state.back();
      _tmp_state.pop_back();
      state = _tmp_state;
    }
    loop->SetState(state);

    // Calculate DEM
    if(loop->parameters.calculate_dem)
    {
      dem->CalculateDEM(i);
    }
    // Save results
    loop->SaveResults(i,time);
    //Update time and counter
    time += tau;
    i++;
  }

  //Set excess number of entries
  int excess = loop->parameters.N - i;
  if(excess<0)
  {
    excess = 0;
  }

  //Print results to file
  loop->PrintToFile(excess);
  if(loop->parameters.calculate_dem)
  {
    dem->PrintToFile(excess);
  }

  //Destroy loop and dem object
  delete loop;
  if(loop->parameters.calculate_dem)
  {
    delete dem;
  }

  return 0;
}
