// ****
// *
// * Header file for miscellaneous functions
// *
// ****

#ifndef MISC_H
#define MISC_H

#include <cstdlib>
#include <vector>

// Find the closest value in a vector to the input x
// Uses a simple search, O(N), so it's best for short arrays
// Returns the index of the vector that's closest to value x
int find_closest(double x, std::vector<double> array, int array_length);

#endif