#ifndef _HYDROGEN_H_
#define _HYDROGEN_H_

#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>
#include <vector>
#include <string>
#include "library.h"
#include "random.h"

class hydrogen {
private:
  double x, y, z; //starting point
  double xtrial, ytrial, ztrial;  // metropolis parameters
  double alpha;
  double step;
  Random rnd;
  int in, tot;

  // possible transition probabilites:
  void uniform_distr();
  void gaussian_distr();

  // distribution probabilities for the first two energy levels:
  double ground_state(double u, double v, double w);
  double exc_state_1(double u, double v, double w);

  // pointers to functions that have double(s) in input, double in output
  // pointer to the correct distribution probability:
  double (hydrogen::*distr_prob)(double x, double y, double z);
  // pointer to the correct transition probability:
  void (hydrogen::*trans_prob)();


public:
  // class constructor (with default values specified ONLY in the .h file!!!)
  // in input:
  // - position coordinates;
  // - step width;
  // - energy level (acceptable: the ground state "100" & first
  //   excited "210")
  // - the transition probability for the Metropolis algorithm,
  //   (acceptable: "uniform" & "gaussian")
  hydrogen(double x=0.,
           double y=0.,
           double z=0.,
           double step=1.,
           std::string energy_level="100",
           std::string tr_prob="uniform");

  ~hydrogen();

  // to reset the position of Metropolis algorithm to sample the prob distr;
  // same inputs as the constructor
  void reset_metropolis(double u,
                   double v,
                   double w,
                   double step=1.,
                   std::string energy_level="100",
                   std::string tr_prob="uniform");

  void metropolis();
  std::vector<double> get_position() const;
  double get_radius() const;
  double acceptance() const;
};


#endif
