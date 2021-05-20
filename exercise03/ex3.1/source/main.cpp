/*
Angela Riva
EXERCISE 03
We compute the value (that is, the price) of EUROPEAN:
- call options;
- put options;
With the hypothesis that the asset price evolution is
S(t)~GBM(mu,sigma^2), geometric brownian motion, and the hypothesis
of a completely efficient market (the evolution of the system depends
only on the situation at t), that is the Markov hypothesis, we calculate
the call and put profit:
- by direct sample of the final price S(T);
- by discretized sample of the GBM path
We estimate the averages and uncertaintes using the blocking method.
We will compare our results with the Black-Scholes analytical solution.
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
  Random rnd;
  rnd.SetRandom("Primes", "seed.in");
  int M = 100000;  // # of estimates
  unsigned int N = 100;  // # of blocks
  int L = int(M/N);  // # of estimates per block
  int Nstep = 100;  // # of discrete time steps
  double S0 = 100.;  // initial price of asset
  double K = 100.;  // strike price
  double T = 1.;  // delivery time
  double r = 0.1; // risk-free interest rate
  double sigma = 0.25;  // volatility
  vector<double> callprice(N);
  vector<double> putprice(N);


  // DIRECT SAMPLING &
  for(unsigned int i=0; i<N; i++) {
    double call=0., put=0.;
    for(int j=0; j<L; j++) {
      double WT = rnd.Gauss(0.,T);
      double ST = S0*exp((r-sigma*sigma*0.5)*T+sigma*WT);
      if(ST>K) call += (ST-K)*exp(-r*T);
      else put += (K-ST)*exp(-r*T);
    }
    callprice[i] = call/double(L);
    putprice[i] = put/double(L);
  }

  vector<double> err = blocking_error(callprice);
  ofstream out("call_direct.dat");
  if(out.fail()){
       std::cerr << "Error opening output file\n";
       return 2;
  }
  for(unsigned int i=0; i<N; i++) {
    out << i << " " << callprice[i] << " " << err[i] << endl;
  }
  out.close();

  out.open("put_direct.dat");
  if(out.fail()){
       std::cerr << "Error opening output file\n";
       return 2;
  }
  err = blocking_error(putprice);
  for(unsigned int i=0; i<N; i++) {
    out << i << " " << putprice[i] << " " << err[i] << endl;
  }
  out.close();

  // DISCRETIZED SAMPLING
  for(unsigned int i=0; i<N; i++) {
    double call=0., put=0.;
    for(int j=0; j<L; j++) {
      double delta = T/double(Nstep);
      double ST = S0;
      for(int k=0; k<Nstep; k++) {
        double Z = rnd.Gauss(0.,1.);
        ST = ST*exp((r-sigma*sigma*0.5)*delta+sigma*Z*sqrt(delta));
      }
      if(ST>K) call += (ST-K)*exp(-r*T);
      else put += (K-ST)*exp(-r*T);
    }
    callprice[i] = call/double(L);
    putprice[i] = put/double(L);
  }

  err = blocking_error(callprice);
  out.open("call_discretized.dat");
  if(out.fail()){
       std::cerr << "Error opening output file\n";
       return 2;
  }
  for(unsigned int i=0; i<N; i++) {
    out << i << " " << callprice[i] << " " << err[i] << endl;
  }
  out.close();

  err = blocking_error(putprice);
  out.open("put_discretized.dat");
  if(out.fail()){
       std::cerr << "Error opening output file\n";
       return 2;
  }
  for(unsigned int i=0; i<N; i++) {
    out << i << " " << putprice[i] << " " << err[i] << endl;
  }
  out.close();

  rnd.SaveSeed();
  return 0;
}
