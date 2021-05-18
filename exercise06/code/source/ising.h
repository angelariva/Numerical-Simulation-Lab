#ifndef __ISING__
#define __ISING__

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <map>
#include <vector>
#include <iomanip>
#include "random.h"

class Ising1D {
private:
  Random rnd;
  // keys for the map corresponding to properties:
  // - energy
  // - capacity
  // - magnetization
  // - susceptibility
  std::vector<std::string> props;
  int seed[4];
  std::map<std::string, double> walker;
  // for blocking method averages:
  std::map<std::string, double> block_average;
  std::map<std::string, double> glob_average;
  std::map<std::string, double> glob_average2;

  double blk_norm, accepted, attempted;
  double stima_u, stima_c, stima_m, stima_x;
  double err_u, err_c, err_m, err_x;

  // configuration
  std::vector<double> s;
  unsigned int n_spin;
  // thermodynamical state
  double beta,temp,J,h;
  // simulation
  int nstep, nblk, metro;
  //output files
  std::ofstream Ene, Heat, Mag, Chi;

  void Input();
  void Reset(unsigned int);
  void Accumulate();
  void BlockAverages(unsigned int);
  void ConfFinal();
  void Measure(int step=0);

  // pointer to be set to the correct sampling method,
  // either Gibbs or Metropolis!
  void (Ising1D::*Move)();
  void MetropolisMove();
  void GibbsMove();

  // inline methods
  double Boltzmann(int, int) const;
  int Pbc(int) const;
  double Error(double,double,int) const;

public:
  // The constructor of the model can start from a previous configuration;
  // the default option is that no configuration file is provided: in this
  // case the values of the spins are selected randomly & the system is
  // equilibrated: the alogrithm is moved for some steps, to allow the
  // temperature to reach a value of equilibrium.
  // the constructor makes use of a file "input.dat" where the Simulation
  // parameters are specified.
  Ising1D(std::string old_configuration="");
  ~Ising1D();

  // run the simulation (default is NOT to print instant values)
  void Run(bool instant=false);
};

/************************************************/
// INLINE METHODS

// Algorithm for periodic boundary conditions
inline int Ising1D::Pbc(int i) const {
  if(i >= (int)n_spin) i = i - n_spin;
  else if(i < 0) i = i + n_spin;
  return i;
}

inline double Ising1D::Error(double sum, double sum2, int iblk) const {
  return std::sqrt((sum2/(double)iblk-pow(sum/(double)iblk,2))/(double)(iblk));
}

inline double Ising1D::Boltzmann(int sm, int ip) const {
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}


#endif
