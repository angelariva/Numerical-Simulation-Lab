#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_
#include <iosfwd>
#include <string>
#include <cmath>
#include <vector>

// function that, given a vector where each element is the average of a block,
// provides, in a vector output:
//  - the progressive blocking error of the vector;
// and modifies the vector input to obtain:
//  - the progressive average of the vector.

std::vector<double> blocking_error(std::vector<double>& average);


// uncertainty calculated as Standard Deviation:
double error(double av, double av2, int n);

#endif
