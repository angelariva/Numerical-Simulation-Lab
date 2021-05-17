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
#include "library.h"

class Ising1D {
private:
  Random rnd;
  // properties:
  // - energy
  // - capacity
  // - magnetization
  // - susceptibility
  vector<std::string> props;
  map<std::string, double> walker;
  // for blocking method averages:
  map<std::string, double> block_average;
  map<std::string, double> glob_average;
  map<std::string, double> glob_average2;

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
  void Measure();

  // pointer to be set to the correct sampling method,
  // either Gibbs or Metropolis!
  void (Ising1D::*Move)();
  void MetropolisMove();
  void GibbsMove();

  // inline methods
  double Boltzmann(int, int);
  int Pbc(int);
  double Error(double,double,int);

public:
  // The constructor of the model can start from a previous configuration;
  // the default option is that no configuration file is provided: in this
  // case the values of the spins are selected randomly.
  Ising1D(std::string old_configuration="");
  ~Ising1D();

  // run the simulation
  void Run();
};

/************************************************/
// INLINE METHODS

// Algorithm for periodic boundary conditions
inline int Ising1D::Pbc(int i) const {
  if(i >= (int)nspin) i = i - nspin;
  else if(i < 0) i = i + nspin;
  return i;
}

inline double Ising1D::Error(double sum, double sum2, int iblk) const {
  if(iblk==1) return 0.0;
  else  return std::sqrt((sum2/(double)iblk-pow(sum/(double)iblk,2))/(double)(iblk-1s));
}

inline double Ising1D::Boltzmann(int sm, int ip) const {
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}


#endif
