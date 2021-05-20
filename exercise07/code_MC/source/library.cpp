#include "library.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>


std::vector<double> blocking_error(std::vector<double> & average) {
  std::vector<double> err(average.size(), 0.);
  double sum=0., sum2=0.;
  sum = average[0];
  sum2 = sum*sum;

  for(unsigned int i=1; i<average.size(); i++) {
    sum += average[i];
    sum2 += average[i]*average[i];
    average[i] = sum/double(i+1);
    err[i] = std::sqrt((sum2/(i+1) - average[i]*average[i])/double(i));
  }
  return err;
}


double error(double av, double av2, int n) {
  if (n==0) return 0;
  else return sqrt((av2 - pow(av,2))/double(n));
}
