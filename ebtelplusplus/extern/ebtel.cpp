/*
ebtel++
A code for computing the evolution of dynamically heated, spatially-averaged solar coronal loops.
*/
#include "boost/numeric/odeint.hpp"
#include "loop.h"
#include "dem.h"
#include "observer.h"


py::dict run(py::dict& config_dict)
{
  //Declarations
  int num_steps;
  state_type state;
  py::dict results_dict;
  LOOP loop;
  DEM dem;
  OBSERVER obs;

  loop = new Loop(config_dict);
  if(loop->parameters.calculate_dem)
  {
    dem = new Dem(loop);
  }
  else
  {
    dem = new Dem();
  }
  obs = new Observer(loop,dem);

  // Set initional conditions of the loop
  state = loop->CalculateInitialConditions();
  // Set initial state for loop and dem
  obs->Observe(state, 0.0);

  // Set up Runge-Kutta integrator
  typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > stepper_type;
  auto controlled_stepper = boost::numeric::odeint::make_controlled(loop->parameters.adaptive_solver_error,
                                                                    loop->parameters.adaptive_solver_error,
                                                                    stepper_type());

  // Integrate
  num_steps = 0;
  if(loop->parameters.use_adaptive_solver)
  {
    // Set maximum number of allowed failures
    int max_failures = 1000;
    // Initialize time and timestep
    double tau = loop->parameters.tau;
    double t = loop->parameters.tau;
    double old_tau,old_t;
    // Start integration loop
    while(t<loop->parameters.total_time)
    {
      int fail = 1;
      int num_failures = 0;
      while(fail>0)
      {
        // Throw error if exceeded max number of failures to avoid infinite loop
        if(num_failures>max_failures)
        {
          throw std::runtime_error("Adaptive solver exceeded maximum number of allowed failures.");
        }
        old_tau = tau;
        old_t = t;
        fail = controlled_stepper.try_step(loop->CalculateDerivatives,state,t,tau);
        // Force NaNs to fail
        if(!fail) fail = obs->CheckNan(state,t,tau,old_t,old_tau);
        num_failures++;
      }
      // Enforce limits on timescale due to thermal conduction and heating
      tau = loop->ControlTimeStep(state, t, tau);
      // Save the state
      obs->Observe(state,t);
      num_steps += 1;
    }
  }
  else
  {
    // Constant timestep integration
    num_steps = boost::numeric::odeint::integrate_const(controlled_stepper,
                                                        loop->CalculateDerivatives,
                                                        state,
                                                        loop->parameters.tau,
                                                        loop->parameters.total_time,
                                                        loop->parameters.tau,
                                                        obs->Observe);
    if(obs->CheckNan(state))
    {
        throw std::runtime_error("NaNs were detected in the output. Check the input configuration.");
    }
  }

  //Package up results into a single Python structure
  if(!loop->parameters.use_adaptive_solver)
  {
    num_steps = std::fmin(loop->parameters.N,num_steps);
  }
  results_dict = loop->GetFinalResults(num_steps);
  if(loop->parameters.calculate_dem)
  {
    results_dict = dem->GetFinalResults(results_dict, num_steps);
  }

  //Cleanup
  delete obs;
  delete loop;
  delete dem;

  return results_dict;
}

PYBIND11_MODULE(_core, m) {
  m.def("run",
        &run,
        py::return_value_policy::copy,
        R"pbdoc(
    Python binding to ebtelplusplus
  )pbdoc");
}
