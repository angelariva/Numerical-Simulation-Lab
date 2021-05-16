/*
Angela Riva
EXERCISE 02.1
Compute the value of a definite integral by:
- sampling a uniform distribution in [0, 1];
- sampling a non uniform distributionin [0, 1]: IMPORTANCE SAMPLING.
We will estimate the averages using the blocking method.
*/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <numeric>
#include "random.h"
#include "library.h"

using namespace std;

// the integrand:
double integrand(double x) {
  return M_PI*0.5*cos(M_PI*0.5*x);
}

// the probability distribution:
double distr(double x) {
  return (1-  M_PI*0.5*(x-0.5));
}

int main (int argc, char *argv[]){
   Random rnd;
   rnd.SetRandom("Primes", "seed.in");
   // Declaring variables:
   unsigned int N = 100;       // Number of blocks
   int M = 100000;   // Number of sample points
   int L = int(M/N);  // Number of sample points in each block

   // EXERCISE 1.1: estimate of the integral with UNIFORM SAMPLING
   vector<double> intgr(N);
   vector<double> eval(M);
   for (auto & el : eval) el = integrand(rnd.Rannyu());

   for(unsigned int i=0; i<N; i++)
     intgr[i] = accumulate(eval.begin()+i*L, eval.begin()+(i+1)*L, 0.)/double(L);
   // error estimation (blocking method)
   vector<double> err = blocking_error(intgr);

   // printing results;
   ofstream out("uniform.dat");
   if(out.fail()){
        std::cerr << "Error opening output file\n";
        return 2;
   }
   for(unsigned int i=0; i<N; i++)
    out << (i+1) << " " << intgr[i] << " " << err[i] << endl;
   out.close();

   // EXERCISE 1.2: ESTIMATE of the integral using IMPORTANCE SAMPLING
   // we use the accept-reject method to sample the chosen probability distribution

   fill(intgr.begin(), intgr.end(), 0.);
   double x; // the point randomly drawn for accept-reject
   double max = 1. + M_PI*0.25; // upper bound of probability distr

   for(unsigned int i=0; i<N; i++){
     for(int j=0; j<L; j++) {
       x = rnd.AcceptReject(distr, max);
       intgr[i] += integrand(x)/distr(x);
     }
     intgr[i] /= double(L);
   }
   std::fill(err.begin(), err.end(), 0);
   err = blocking_error(intgr);

   // printing results;
   out.open("importance.dat");
   if(out.fail()){
        std::cerr << "Error opening output file\n";
        return 2;
   }
   for(unsigned int i=0; i<N; i++) {
    out << (i+1) << " " << intgr[i] << " " << err[i] << endl;
   }
   out.close();
   rnd.SaveSeed();
   return 0;
}
