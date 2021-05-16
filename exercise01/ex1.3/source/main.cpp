/*
Angela Riva

01.3) Buffon's experiment
To get an estimate of the value of Pi we simulate Buffon's experimet.
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
  double l=0.5, d=1.; // length of needle, separation of 'tiles'
  double av=0., av2=0., sum=0., sum2=0., err=0.;
  int hit=0; // counters
  // # of throws, # of blocks, # throws in each block
  int M=100000, N=100, L=int(M/N);
  ofstream out("results3.dat");

  for(int i=0; i<N; i++){
    hit=0;
    for(int j=0; j<L; j++){
      double x=rnd.Rannyu(0., d);
      double x1,y1,r;
      do{
        x1 = rnd.Rannyu();
        y1 = rnd.Rannyu();
        r = sqrt(x1*x1 + y1*y1);
      }while(r>=1.);
      if((x-l/2.*x1/r)<=0. || (x+l/2.*x1/r)>=d) hit++;
    }
    av = 2*l*L/(double)(hit*d);
    av2 = pow(av,2);
    sum = (sum*i+av)/(double)(i+1);
    sum2 = (sum2*i+av2)/(double)(i+1);
    if(i==0) err = 0;
    else err = error(sum, sum2, i);

    out << (i+1) << " " << sum << " " << err << endl;
  }

  out.close();
  rnd.SaveSeed();
  return 0;
}
