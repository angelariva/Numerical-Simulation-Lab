#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <map>
#include <vector>
#include <iostream>
#include "random.h"

class MolDynMC {
private:
  Random rnd;
  int seed[4]; // to initialize the random number generator

  // coordinates for molecules
  std::vector<double> x, y, z;
  // keys for the map corresponding to observables:
  // - energy
  // - pressure
  std::vector<std::string> props;
  // for blocking method averages (of properties):
  std::map<std::string, double> walker;
  std::map<std::string, double> block_average;
  std::map<std::string, double> glob_average;
  std::map<std::string, double> glob_average2;

  // to compute g(r) (you need histograms):
  // for blocking method averages:
  std::vector<double> histo_walker;
  std::vector<double> histo_block_average;
  std::vector<double> histo_glob_average;
  std::vector<double> histo_glob_average2;

  // parameters:
  double blk_norm, accepted, attempted;
  double box, temp,  rho, rcut, delta;
  unsigned int npart, nblk, nstep;
  double beta, vol, vtail, ptail;
  const unsigned int nbins = 100;
  double bin_size;

  std::ofstream Gerr, Gave, Epot, Pres;
  std::ofstream ist_pot;
  std::ofstream ist_pres;

  bool instant;
  // METHODS
  void Reset(unsigned int); // to reset the block average
  void Accumulate();  // get block average
  void Averages(unsigned int);  // blocking averages
  void Move(); // one MC step
  void ConfFinal();
  void ConfXYZ(unsigned int);
  void Measure(unsigned int istep=0);
  double Boltzmann(double, double, double, unsigned int);
  double Pbc(double) const; // periodic boundary conditions
  double Error(double,double,unsigned int) const;   // uncertainties

public:
  // CONSTRUCTOR!!
  // - 'simParameters' is the filename for input parameters
  // - "initial_configuration" is the filename of a file containing
  //   a configuration valid for the molecular model defined in
  //   "input.dat".
  // - istantaneous: if true the program will print istantaneous value of
  //                 p and U, and not just the mean with data blocking.
  MolDynMC(std::string simParameters, std::string initial_configuration,
              bool instantaneous=false);

  //Run the simulation
  void Run();
  ~MolDynMC();

};


//***************************************************************//
// inline METHODS************************************************//

inline double MolDynMC::Pbc(double r) const {
  return r - box * rint(r/box);
}

inline double MolDynMC::Error(double sum, double sum2, unsigned int iblk) const {
  return std::sqrt((sum2/(double)iblk - std::pow(sum/(double)iblk,2))/(double)iblk);
}


#endif
