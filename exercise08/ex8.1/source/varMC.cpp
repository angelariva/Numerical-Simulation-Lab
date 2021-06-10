#include "varMC.h"

using namespace std;

Particle::Particle(double mu, double sigma, double x_start, double x_jump) :
    mu(mu),
    sigma(sigma),
    x(x_start),
    x_jump(x_jump) {
      rnd.SetRandom("Primes", "seed.in");
      in = 0;
      tot = 0;
      mu_opt = NAN;
      sigma_opt = NAN;
}

void Particle::Reset(double m, double s, double x_start, double x_j) {
  mu = m;
  sigma = s;
  x = x_start;
  x_jump = x_j;
  in = 0;
  tot = 0;
}

double Particle::Distribution(double a) {
  double b;
  b = exp(-(a-mu)*(a-mu)/(sigma*sigma*2.)) + exp(-(a+mu)*(a+mu)/(sigma*sigma*2.));
  return b*b;
}

double Particle::Sample() {
  Metropolis();
  return x;
}

void Particle::Metropolis() {
    x_trial = (rnd.Rannyu()-0.5)*2.*x_jump + x;
    if(rnd.Rannyu()<min(1., Distribution(x_trial)/Distribution(x))) {
        in++;
        x = x_trial;
    }
    tot++;
}

double Particle::EnergyLoc(double a) {
  double s, T, V;
  s = sigma*sigma;
  T = -0.5/s*((a*a + mu*mu)/s - 1. - (2.*a*mu/s)*tanh(a*mu/s));
  V = (a*a - 2.5)*a*a;
  return T + V;
}

double Particle::Acceptance() {
  return double(in)/tot;
}

void Particle::Results(double & mu, double & sigma) {
  mu = mu_opt;
  sigma = sigma_opt;
}
