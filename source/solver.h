/*
solver.h
Class definition for 1D ODE solvers
*/

#ifndef SOLVER_H
#define SOLVER_H

#include "helper.h"
#include "loop.h"

// Solver object
//
// Class for holding various solvers. These include an Euler solver and a Fourth-order
// Runge-Kutta solver in addition to an adaptive timestep routine to be used in conjunction
// with the Runge-Kutta solver.
//
class Solver {
private:

  LOOP loop;

public:

  // Constructor for <Solver> class.
  // @loop_object and instance of <Loop>.
  //
  // Set up the solver class and declare <loop_object> as a member of <Solver>.
  //
  Solver(LOOP loop_object);

  // Destructor
  ~Solver();

  // Euler solver
  // @state current state of the loop
  // @time current time (in s)
  // @tau timestep (in s)
  //
  // Evolve <state> in time using a forward Euler method by a timestep <tau>.
  //
  // @return the evolved loop state
  //
  std::vector<double> EulerSolver(std::vector<double> state, double time, double tau);

  // Fourth-order Runge-Kutta solver
  // @state current state of the loop
  // @time current time (in s)
  // @tau timestep (in s)
  //
  // Evolve <state> in time using a basic fourth-order Runge-Kutta method.
  //
  // @return the evolved loop state
  //
  std::vector<double> RK4Solver(std::vector<double> state, double time, double tau);

  // Adaptive time-stepper routine
  // @state current state of the loop
  // @time current time (in s)
  // @tau current timestep (in s)
  //
  // Evolve <state> in time using a fourth-order Runge-Kutta method, but using a step-doubling
  // method where the truncation error between a step using <tau> and <tau>/2 must be below a
  // given error tolerance as set in the configuration file. This routine is similar to that outlined
  // in [Garcia (1994)](http://www.algarcia.org/nummeth/nummeth.html).
  //
  // @return the evolved loop state
  //
  std::vector<double> RKA4Solver(std::vector<double> state, double time, double tau);
};
// Pointer to the <Solver> class
typedef Solver* SOLVER;

#endif
