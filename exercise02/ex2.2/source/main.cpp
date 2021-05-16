/*
Angela Riva
EXERCISE 02.2
Simulate a Random Walk starting from the origin:
- on a 3D lattice (discrete);
- on a random direction in the solid angle 4pi (continum);
Calculate the mean positions, their uncertaintes as functions of the STEPS.
// N.B. THE STEP OF THE LATTICE IS a=1!
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
#include "walk.h"

using namespace std;


int main (int argc, char *argv[]){
   Random rnd;
   rnd.SetRandom("Primes", "seed.in");
   unsigned int N = 100;       // Number of blocks
   int M = 10000;   // Number of random walks
   int L = int(M/N);  // Number of random walks in each block
   unsigned int Nstep = 100; // Number of steps in each random walk

   // EXERCISE 2.1: simulation of random walk on discrete 3D lattice
   ofstream out("discrete.dat");
   ofstream out1("discrete_blocking.dat");
   if(out.fail()){
        std::cerr << "Error opening output file\n";
        return 2;
   }
   // loop on the # of steps
   for(unsigned int i=0; i<Nstep; i++) {
     vector<double> r2Discrete(N);  // it will contain the mean r^2 for each block
     // loop on the N of blocks
     for(unsigned int j=0; j<N; j++) {
       // a loop on the L random walks
       for(int k=0; k<L; k++){
         // Memory is dinamically allocated here to exploit c++ polymorphism:
         // we want to use the inherited method 'walk'
         RandomWalk* w = new DiscreteW();
         for(unsigned int l=0; l<i; l++){
           w->walk(&rnd); // w makes i steps...
         }
         r2Discrete[j] += pow(w->get_dist(), 2);
         // we deallocate the memory!
         delete w;
       }
       r2Discrete[j] /= double(L);
     }
     vector<double> err = blocking_error(r2Discrete);
     out << i << " " << sqrt(r2Discrete[N-1]) << " " << err[N-1] << endl;
     if(i==Nstep-1) {
       for(unsigned int k=0; k<N; k++)
        out1 << k << " " << r2Discrete[k] << " " << err[k] << endl;
     }
   }

   out.close();
   out1.close();

   // EXERCISE 2.2: simulation of random walk on random 3D direction
  out.open("continuous.dat");
  out1.open("continuous_blocking.dat");

   if(out.fail()){
        std::cerr << "Error opening output file\n";
        return 2;
   }
   // loop on the # of steps
   for(unsigned int i=0; i<Nstep; i++) {
     vector<double> r2Continuous(N);  // it will contain the mean r^2 for each block
     // loop on the N of blocks
     for(unsigned int j=0; j<N; j++) {
       // a loop on the L random walks
       for(int k=0; k<L; k++){
         // Memory is dinamically allocated here to exploit c++ polymorphism:
         // we want to use the inherited method 'walk'
         RandomWalk* w = new ContinousW();
         for(unsigned int l=0; l<i; l++){
           w->walk(&rnd); // w makes i steps...
         }
         r2Continuous[j] += pow(w->get_dist(), 2);
         // we deallocate the memory!
         delete w;
       }
       r2Continuous[j] /= double(L);
     }
     vector<double> err = blocking_error(r2Continuous);
     out << i << " " << sqrt(r2Continuous[N-1]) << " " << err[N-1] << endl;
     if(i==Nstep-1) {
       for(unsigned int k=0; k<N; k++)
        out1 << k << " " << r2Continuous[k] << " " << err[k] << endl;
     }
   }

   out.close();
   out1.close();

   rnd.SaveSeed();
   return 0;
}
