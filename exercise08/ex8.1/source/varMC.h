#ifndef _VAR_MC_H_
#define _VAR_MC_H_

#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "random.h"

class Particle {
private:
  // variational parameters to set the wavefunction:
  double mu, sigma;
  double mu_opt, sigma_opt;
  // Metropolis parameters:
  double x, x_trial, x_jump;
  int in, tot;
  Random rnd;

public:
  Particle(double mu, double sigma, double x_start=0., double x_jump=1.);
  ~Particle() { };

  void Reset(double mu, double sigma, double x_start=0., double x_jump=1.);
  void Metropolis();                  // Metropolis algorithm computed with
                                      // UNIFORM transition probability, samples
                                      // psi^2, the density distribution function of the particle
  double Sample();                // Executes Metropolis once, returns the x value
  double Distribution(double a);  // probability density give by the wavefunction
  double Acceptance();            // acceptance rate of Metropolis Algorithm
  double EnergyLoc(double a);     // Local Energy: (T+V)psi/psi, hbar = m = 1.
                                  // T = -d^2/dx^2
  void Results(double & mu, double & sigma);
};



#endif //_VAR_MC_H_
