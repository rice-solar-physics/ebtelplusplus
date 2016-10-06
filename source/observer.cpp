/* observer.cpp
Function definitions for Observer methods
*/

#include "observer.h"

int Observer::i;
LOOP Observer::loop;
DEM Observer::dem;

Observer::Observer(LOOP loop_object,DEM dem_object)
{
  // Initialize counter
  i = 1;
  // Set needed objects
  loop = loop_object;
  dem = dem_object;
}

Observer::~Observer(void)
{
  // Destroy things here
}

void Observer::Observe(const state_type &state, const double time)
{
  // Store state
  loop->SetState(state);
  // Calculate DEM
  if(loop->parameters.calculate_dem)
  {
    dem->CalculateDEM(i);
  }
  // Save results
  loop->SaveResults(i,time);
  // Increment counter
  i++;
}
