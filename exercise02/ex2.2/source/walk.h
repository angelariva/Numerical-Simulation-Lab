#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "random.h"


#ifndef __Walk__
#define __Walk__

class RandomWalk {
protected:
  int Nstep;
  double dist;  // distance from the origin
  virtual void set_dist() = 0;
public:
  RandomWalk();
  virtual ~RandomWalk();
  double get_dist();
  int get_step();

  // method that any class inherits and implpements: how to make a step
  virtual void walk(Random*) = 0;
};

class DiscreteW : public RandomWalk {
private:
  int pos[3];   // the position is given by three coordinates
  virtual void set_dist();
public:
  int get_x();
  int get_y();
  int get_z();
  DiscreteW();
  ~DiscreteW();
  virtual void walk(Random* rnd);
};

class ContinousW : public RandomWalk {
private:
  double pos[3];
  virtual void set_dist();
public:
  ContinousW();
  ~ContinousW();
  virtual void walk(Random* rnd);
};

#endif // __Walk__
