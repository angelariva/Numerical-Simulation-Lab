/*
Angela Riva
2.1) Implementation of NEW METHODS for the Pseudo-Random Numbers Generator:
- Exponential distribution
- Cauchy-Lorentz distribution
2.2) TEST the validity of the Central Limit Theorem
by exctracting and averaging for 10^4 times N=1,2,10,100 random variables.
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

int main (int argc, char *argv[]){

  // Declaring variables:
  Random rnd;
  rnd.SetRandom("Primes", "seed.in");
  int M = 10000;
  double mu=0., gamm=1., lam=1.;
  vector<int> N = {1,2,10,100};
  ofstream out1("results21.dat");
  ofstream out2("results22.dat");
  ofstream out3("results23.dat");


  // A loop to repeat the experiment M times:
  for(int k=0; k<4; k++){
    for(int i=0; i<M; i++){
      double gaus=0., exp=0., lor=0.;
      // A loop to get the average of N random variables:
      for(int j=0; j<N[k]; j++){
          gaus += rnd.Rannyu();
          exp += rnd.Exp(lam);
          lor += rnd.Lorentz(gamm, mu);
      }
      gaus /= double(N[k]);
      exp /= double(N[k]);
      lor /= double(N[k]);

      out1 << gaus << endl;
      out2 << exp << endl;
      out3 << lor << endl;
    }
  }

  out1.close();
  out2.close();
  out3.close();

  return 0;
}
