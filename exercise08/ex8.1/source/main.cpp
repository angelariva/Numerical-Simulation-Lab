/*
Angela Riva
EXERCISE 08

08.1 & 08.2)
      Find the Ground State Energy of a single 1D particle with:
      H = p^2/2m + x^4 - 2.5*x^2,
      Use a variational MC approach starting from a trial wavefunction
      made up by the "superposition" of two gaussians centered in mu, -mu,
      whith sigma as variance.
      sigma and mu are the variational parameters, and they will be chosen
      as the minimizers of the GS energy calculated with data blocking.
      The integral is calculated using the Metropolis algorithm, with a
      uniform transition probability.
*/



#include <vector>
#include <string>
#include "random.h"
#include "varMC.h"
#include "library.h"

using namespace std;

int main(int argc, char *argv[]) {

  std::string usage =
  "Usage: " + std::string(argv[0]) + " [mu] [sigma] [M]\n" +
  "This program computes the GS energy of a 1D particle using Metropolis\n"+
  "the hamiltonian is: \n"+
  "\t H = p^2/2m + x^4 - 2.5*x^2 \n"+
  "the trial wavefunction is: \n"+
  "\t psi(x) ~ e^{(x-mu)^2/2sigma^2} + e^{(x+mu)^2/2sigma^2} \n"+
  "The program accepts in input sigma and mu as parameters (MANDATORY).\n"+
  "[M] is the optional parameter that specifies the number of points to sample.";

  if(argc!=3 and argc!=4) cout << usage << endl;
  else {
    double mu;
    double sigma;

    mu = atof(argv[1]);
    sigma = atof(argv[2]);


    cout << "mu = " << mu << ", " << argv[1] << ", sigma = " << sigma << ", " << argv[2]  << endl;
    Particle p(mu, sigma, 0., 2.7);
    ofstream outen("energy.dat");
    ofstream outpsi("psi2.dat");
    int M;  // # of throws
    int N = 100;    // # of blocks
    if (argc==3) M = 100000;
    else M = atof(argv[3]);
    int L = int(M/N);    // # of throws in each block

    std::vector<double> psi2;

    double a, sum=0., sum2=0., err=0.;

    for(int i=0; i<N; i++){
      double en=0., en2=0.;
      for(int j=0; j<L; j++){
        a = p.Sample();
        en += p.EnergyLoc(a);
        psi2.push_back(a);
      }
      en = en/double(L);
      en2 = en*en;
      sum = (sum*i + en)/(double(i+1.));
    	sum2 = (sum2*i + en2)/(double(i+1.));
      err = error(sum, sum2, i);
      outen << i+1 << " " << sum << " " << err << endl;
    }
    cout << "Acceptance rate of Metropolis alg: " << p.Acceptance() << endl;

    for(int i=0; i<L; i++) {
      outpsi << psi2[i] << endl;
    }
  }

  return 0;
}
