#include "walk.h"

using namespace std;

RandomWalk :: RandomWalk() {Nstep=0, dist=0;}
RandomWalk :: ~RandomWalk() { }
double RandomWalk :: get_dist() { return dist;}
int RandomWalk :: get_step() { return Nstep;}


void DiscreteW :: set_dist() {
  double d2=0;
  for(auto& el : pos) d2 += el*el;
  dist = sqrt(d2);
}

int DiscreteW :: get_x() { return pos[0]; }
int DiscreteW :: get_y() { return pos[1]; }
int DiscreteW :: get_z() { return pos[2]; }

DiscreteW :: DiscreteW() { for(auto& el : pos) el = 0; }
DiscreteW :: ~DiscreteW() { }

void DiscreteW :: walk(Random* rnd) {
  int a = (int)rnd->Rannyu(1,7);
  if(a%2==0) pos[a/2-1] += 1;
  else pos[(a+1)/2-1] += -1;
  set_dist();
  Nstep++;
}

void ContinousW :: set_dist() {
  double d2=0;
  for(auto& el : pos) d2 += el*el;
  dist = sqrt(d2);
}
ContinousW :: ContinousW() { for(auto& el : pos) el = 0.; }
ContinousW :: ~ContinousW() { }
void ContinousW :: walk(Random* rnd) {
  double phi, theta;
  rnd->SolidAngle(theta, phi);
  pos[0] += sin(theta)*cos(phi);
  pos[1] += sin(theta)*sin(phi);
  pos[2] += cos(theta);
  set_dist();
  Nstep++;
}
