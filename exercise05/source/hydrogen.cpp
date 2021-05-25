#include "hydrogen.h"

hydrogen::hydrogen(double u,
         double v,
         double w,
         double stepp,
         std::string energy_level,
         std::string tr_prob) :
         x(u), y(v), z(w), step(stepp) {

    if (energy_level=="100") distr_prob = &hydrogen::ground_state;
    if (energy_level=="210") distr_prob = &hydrogen::exc_state_1;

    if (tr_prob=="uniform") trans_prob = &hydrogen::uniform_distr;
    if (tr_prob=="gaussian") trans_prob = &hydrogen::gaussian_distr;

    rnd.SetRandom("Primes", "seed.in");
    in=0;
    tot=0;
}

hydrogen::~hydrogen() { }

void hydrogen::reset_metropolis(double u,
         double v,
         double w,
         double stepp,
         std::string energy_level,
         std::string tr_prob) {

    x = u;
    y = v;
    x = w;
    step = stepp;
    if (energy_level=="100") distr_prob = &hydrogen::ground_state;
    if (energy_level=="210") distr_prob = &hydrogen::exc_state_1;

    if (tr_prob=="uniform") trans_prob = &hydrogen::uniform_distr;
    if (tr_prob=="gaussian") trans_prob = &hydrogen::gaussian_distr;

    in=0;
    tot=0;
    xtrial=0.;
    ytrial=0.;
    ztrial=0.;
    alpha=0.;
    rnd.SetRandom("Primes", "seed.in");

}

void hydrogen::metropolis() {
  (*this.*trans_prob)();
  alpha = std::min(1.,
    ((*this.*distr_prob)(xtrial, ytrial, ztrial))/((*this.*distr_prob)(x,y,z)));

  if (rnd.Rannyu() <= alpha) { //valid even if alpha = 1
    x = xtrial;
    y = ytrial;
    z = ztrial;
    in++;
  }
  tot++;
}

std::vector<double> hydrogen::get_position() const {
  std::vector<double> pos={x,y,z};
  return pos;
}

double hydrogen::get_radius() const {
  return std::sqrt(x*x+y*y+z*z);
}

double hydrogen::acceptance() const {
  return double(in)/tot;
}

double hydrogen::ground_state(double u, double v, double w) {
  return std::exp(-2.*std::sqrt(u*u + v*v + w*w))/M_PI;
}

double hydrogen::exc_state_1(double u, double v, double w) {
  double r = u*u + v*v + w*w;
  return std::exp(-std::sqrt(r))*r*r/(32.*M_PI);
}

void hydrogen::uniform_distr() {
  xtrial = x + rnd.Rannyu(-step, step);
  ytrial = y + rnd.Rannyu(-step, step);
  ztrial = z + rnd.Rannyu(-step, step);
}

void hydrogen::gaussian_distr() {
  xtrial = rnd.Gauss(x, step/2.);
  ytrial = rnd.Gauss(y, step/2.);
  ztrial = rnd.Gauss(z, step/2.);
}
