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

/*
void Particle::Optimize(double mu_start,
                        double mu_var,
                        double sigma_start,
                        double sigma_var,
                        double x_start,
                        double x_j,
                        int n_samples,
                        int N_tot,
                        int iter)
                        {
  double mu_min, mu_delta, mu_trial, mu_tmp=0.;
  double sigma_min, sigma_delta, sigma_trial;
  double integral = 1E10;
  double integral_old = 0.;

  mu_min = mu_start;
  sigma_min = sigma_start;

  for(int i=0; i<iter; i++){
    sigma_delta = 2.*sigma_var/double(n_samples*pow(10., i));
    sigma_trial = sigma_min - sigma_var/pow(10., i);
    mu_delta = 2.*mu_var/double(n_samples*pow(10., i));
    mu_trial = mu_min - mu_var/pow(10., i);
    mu_tmp = mu_trial;

    for(int k=0; k<n_samples; k++) {
            for(int z=0; z<n_samples; z++) {
                Reset(mu_trial, sigma_trial, x_start, x_j);
                for(int j=0; j<N_tot; j++) integral += Sample();
                integral /= N_tot;
                if(integral < integral_old) {
                    mu_min = mu_trial;
                    sigma_min = sigma_trial;
                }
                integral_old = integral;
                integral = 0.;
                mu_trial += mu_delta;
            }
            mu_trial = mu_tmp;
            sigma_trial += sigma_delta;
        }

  }
  mu_opt = mu_min;
  sigma_opt = sigma_min;
}
*/

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
