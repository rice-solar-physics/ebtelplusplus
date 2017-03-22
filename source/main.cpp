/*
ebtel++
A code for computing the evolution of dynamically heated, spatially-averaged solar coronal loops.
*/

#include <time.h>
#include "boost/program_options.hpp"
#include "boost/numeric/odeint.hpp"
#include "loop.h"
#include "dem.h"
#include "observer.h"

int main(int argc, char *argv[])
{
  //Declarations
  int num_steps;
  state_type state;
  char rad_config[256],ebtel_config[256];
  LOOP loop;
  DEM dem;
  OBSERVER obs;

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

  // Create loop object
  loop = new Loop(ebtel_config,rad_config);
  // Create DEM object
  if(loop->parameters.calculate_dem)
  {
    dem = new Dem(loop);
  }
  else
  {
    dem = new Dem();
  }
  // Configure observer
  obs = new Observer(loop,dem);

  // Set initional conditions of the loop
  state = loop->CalculateInitialConditions();
  // Set initial state for loop and dem
  obs->Observe(state, 0.0);

  // Set up Runge-Kutta integrator
  typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > stepper_type;
  auto controlled_stepper = boost::numeric::odeint::make_controlled(loop->parameters.adaptive_solver_error, loop->parameters.adaptive_solver_error, stepper_type());
  // Integrate
  num_steps = 0;
  if(loop->parameters.use_adaptive_solver)
  {
    // Initialize time and timestep
    double tau = loop->parameters.tau;
    double t = loop->parameters.tau;
    // Start integration loop
    while(t<loop->parameters.total_time)
    {
      int fail = 1;
      while(fail)
      {
        fail = controlled_stepper.try_step(loop->CalculateDerivs,state,t,tau);
        if(!fail)
        {
          fail = obs->CheckNan(state,t,tau);
        }
      }
      // Enforce thermal conduction timescale limit
      double tau_tc = 4e-10*state[2]*pow(loop->parameters.loop_length,2)*pow(std::fmax(state[3],state[4]),-2.5);
      tau = std::fmin(tau,0.5*tau_tc);
      // Enforce limit set by duration of heating events
      tau = std::fmin(tau,loop->GetMaxAllowedTimestep());
      obs->Observe(state,t);
      num_steps += 1;
    }
  }
  else
  {
    num_steps = boost::numeric::odeint::integrate_const( controlled_stepper, loop->CalculateDerivs, state, loop->parameters.tau, loop->parameters.total_time, loop->parameters.tau, obs->Observe);
  }

  //Print results to file
  if(!loop->parameters.use_adaptive_solver)
  {
    num_steps = std::fmin(loop->parameters.N,num_steps);
  }
  loop->PrintToFile(num_steps);
  if(loop->parameters.calculate_dem)
  {
    dem->PrintToFile(num_steps);
  }

  //Cleanup
  delete obs;
  delete loop;
  delete dem;

  return 0;
}
