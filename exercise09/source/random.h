/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SetRandom(std::string, std::string);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  int RandInt(int min, int max);
  double Gauss(double mean, double sigma);
  double Exp(double lambda);
  double Lorentz(double gamm, double mu);

  // Angle: generates a uniform distribution given angle and radius!
  double Angle();

  // SolidAngle: void function. the arguments are passed by reference
  // the function modyfies them in order to generate a uniform distribution
  // of points on a sphere centered in the origin (radius is not important)
  // the distribution is obtained with inversion method
  void SolidAngle(double& theta, double& phi);

  // AcceptReject: method to generate random numbers distributed as f in [0,1).
  // f is a pointer to a function, specifies the desired analytical form ouf
  // the distribution;
  // max is the maximum value of f in the interval considered
  // the output is a random variable x distributed as f.
  double AcceptReject(double (*f)(double), double max=1);
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
