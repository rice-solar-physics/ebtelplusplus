// ****
// *
// * Header file for miscellaneous functions
// *
// ****

#ifndef MISC_H
#define MISC_H

#include <cstdlib>
#include <vector>

// Find the closest value in an array (or vector) to the input x
// Uses a simple search, O(N), so it's best for short arrays
// Returns the index of the array (vector) that's closest to value x
int find_closest(double x, double array[], int array_length);
int find_closest(double x, std::vector<double> array, int array_length);

#endif