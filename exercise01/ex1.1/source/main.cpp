/*
Angela Riva

TEST of the Pseudo-Random Numbers Generator:
Estimate of:
1.1) <r> and of its uncertainty;
1.2) <sigma^2> and of its uncertainty;
both with blocking method (as functions of # of blocks).
1.3) chi^2 test
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
   // Declaring variables:
   unsigned int N = 100;       // Number of blocks
   int M = 1000000;   // Number of throws
   int L = int(M/N);  // Number of throws in each block
   vector<double>  distr(M);
   for(auto & elem : distr) elem = rnd.Rannyu();

   /************************************************************/
   // EXERCISE 1.1: <r> and uncertainty;
   vector<double> av(N);

   // getting the averages of the N blocks:
   for(unsigned int i=0; i<N; i++) {
     // accumulate(first_vector_elem_to_sum, last_vector_elem_to_sum, initial_sum_value)
     av[i] = accumulate(distr.begin()+i*L, distr.begin()+(i+1)*L, 0.)/double(L);
   }

   // get the average and the error as a function of the # of the blocks:
   vector<double> error = blocking_error(av);

   // printing results;
   ofstream out("results11.dat");
   if(out.fail()){
        std::cerr << "Error opening output file\n";
        return 2;
   }
   for(unsigned int i=0; i<N; i++)
    out << (i+1) << " " << av[i] << " " << error[i] << endl;

   out.close();


   /************************************************************/
   // EXERCISE 1.2: <sigma^2> and uncertainty;
   std::fill(av.begin(), av.end(), 0);
   std::fill(error.begin(), error.end(), 0);

   // We proceed as in the first part of the exercise;
   // this time we will compute <(r-0.5)^2> (*):
   for(unsigned int i=0; i<N; i++) {
     // using a lambda function in accumulate to avoid double loop;
     // the lambda computes (*):
     av[i] = accumulate(distr.begin()+i*L, distr.begin()+(i+1)*L, 0.,
        [](double lhs, double rhs){return lhs+(rhs-0.5)*(rhs-0.5);})/double(L);
   }
   error = blocking_error(av);

   // printing results;
   out.open("results12.dat");
   if(out.fail()){
     std::cerr << "Error opening output file\n";
     return 2;
   }
   for(unsigned int i=0; i<N; i++)
    out << (i+1) << " " << av[i] << " " << error[i] << endl;
   out.close();

   /************************************************************/
   // EXERCISE 1.3: chi^2 test;

   N = 100; // # of subintervals
   M = 10000;   // # of throws
   L = int(M/N); // expected number of throws in each subinterval

   distr.resize(N*M);
   for(auto & elem : distr) elem = rnd.Rannyu();
   out.open("results13.dat");

   vector<int> events(N); // vector containing the # of events in each interval
   double chisq=0;
   // We compute chi^2 100 times:
   for(int j=0; j<100; j++) {
     std::fill(events.begin(), events.end(), 0);
     for(int i=0; i<M; i++) events[int(distr[j*M+i]*N)]++;
     for(unsigned int i=0; i<N; i++) chisq += pow(events[i] - L, 2)/L;
     out << chisq << endl;
     chisq = 0;
   }

   out.close();

   rnd.SaveSeed();
   return 0;
}
