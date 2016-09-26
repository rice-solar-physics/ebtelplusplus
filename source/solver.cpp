/*
solver.cpp
Methods and attributes for Solver class
*/

#include "solver.h"

Solver::Solver(LOOP loop_object)
{
  loop = loop_object;
}

Solver::~Solver(void)
{
  // Free some stuff here if needed
}

std::vector<double> Solver::EulerSolver(std::vector<double> state,double time,double tau)
{
  std::vector<double> new_state(state.size());
  std::vector<double> derivs = loop->CalculateDerivs(state, time);

  for(int i = 0;i<state.size();i++)
  {
    new_state[i] = state[i] + derivs[i]*tau;
  }

  return new_state;
}

std::vector<double> Solver::RK4Solver(std::vector<double> state, double time, double tau)
{
  int i;
  std::vector<double> new_state(state.size());
  std::vector<double> _tmp_state(state.size());
  std::vector<double> f1,f2,f3,f4;

  f1 = loop->CalculateDerivs(state,time);
  for(i=0;i<f1.size();i++)
  {
    _tmp_state[i] = state[i] + f1[i]*tau/2.0;
  }

  f2 = loop->CalculateDerivs(_tmp_state,time+tau/2.0);
  for(i=0;i<f2.size();i++)
  {
    _tmp_state[i] = state[i] + f2[i]*tau/2.0;
  }

  f3 = loop->CalculateDerivs(_tmp_state,time+tau/2.0);
  for(i=0;i<f3.size();i++)
  {
    _tmp_state[i] = state[i] + f3[i]*tau;
  }

  f4 = loop->CalculateDerivs(_tmp_state,time+tau);

  for(i=0;i<f4.size();i++)
  {
    new_state[i] = state[i] + (f1[i] + f4[i] + 2.0*(f2[i] + f3[i]))*tau/6.0;
  }

  return new_state;
}

std::vector<double> Solver::RKA4Solver(std::vector<double> state, double time, double tau)
{
  std::vector<double> small_step;
  std::vector<double> big_step;
  std::vector<double> result;
  std::vector<double> _tmp_error_ratio(state.size());
  double scale,diff,old_tau;
  double error_ratio;

  int i = 0;
  int max_try = 100;
  double safety_1 = 0.9;
  double safety_2 = 1.1;
  double safety_3 = 4.0;
  double epsilon = 1.0e-16;

  for(i=0;i < max_try;i++)
  {
    // Two small steps
    small_step = RK4Solver(state,time,0.5*tau);
    small_step = RK4Solver(small_step,time+0.5*tau,0.5*tau);
    // Big step
    big_step = RK4Solver(state,time,tau);
    // Calculate error ratio
    for(int j=0;j<_tmp_error_ratio.size();j++)
    {
      scale = loop->parameters.rka_error*(std::abs(small_step[j]) + std::abs(big_step[j]))/2.0;
      diff = small_step[j] - big_step[j];
      _tmp_error_ratio[j] = std::abs(diff)/(scale+epsilon);
    }
    error_ratio = *std::max_element(_tmp_error_ratio.begin(),_tmp_error_ratio.end());
    // Estimate new value of tau
    old_tau = tau;
    tau = std::fmax(safety_1*tau*std::pow(error_ratio,-1.0/5.0),tau/safety_2);
    if(error_ratio<1)
    {
      break;
    }
  }

  if(i==max_try)
  {
    std::cout << "Warning! Adaptive solver did not converge to best step size." << std::endl;
  }

  //TODO: add thermal conduction timestep limitation
  // Limit timestep increase
  tau = std::fmin(tau,safety_3*old_tau);
  // Add the timestep to the returned state
  small_step.push_back(tau);

  return small_step;
}
