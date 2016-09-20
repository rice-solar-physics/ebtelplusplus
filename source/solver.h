/*
solver.h
Class definition for 1D ODE solvers
*/

#ifndef SOLVER_H
#define SOLVER_H

#include "helper.h"
#include "loop.h"

class Solver {
private:
  LOOP loop;
public:
  // Constructor for <Solver> class.
  // @loop_object and instance of <Loop>.
  //
  Solver(LOOP loop_object);

  // Destructor
  ~Solver();

  // Euler solver
  //
  std::vector<double> EulerSolver(std::vector<double> state, double time, double tau);

  // Fourth-order Runge-Kutta solver
  //
  std::vector<double> RK4Solver(std::vector<double> state, double time, double tau);

  // Adaptive time-stepper for RK4 solver
  //
  std::vector<double> RKA4Solver(std::vector<double> state, double time, double tau);
};
typedef Solver* SOLVER;

#endif
