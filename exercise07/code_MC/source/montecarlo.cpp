#include "montecarlo.h"

MolDynMC::MolDynMC(std::string simParameters, std::string initial_configuration, bool instantaneous) : instant(instantaneous) {

  props = {"energy", "pressure"};

  std::ifstream ReadInput(simParameters);
  if(ReadInput.fail()) {
    std::cerr << "No valid " << simParameters << " file found. Aborting.\n";
    exit(1);
  }

  std::ifstream ReadConf(initial_configuration);
  if(ReadConf.fail()) {
    std::cerr << "Confifuration file " << initial_configuration
    << " not found. Aborting.\n";
    exit(1);
  }

  std::cout << "Classic Lennard-Jones fluid        " << std::endl;
  std::cout << "Monte Carlo simulation             " << std::endl << std::endl;
  std::cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << std::endl << std::endl;
  std::cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << std::endl << std::endl;
  std::cout << "The program uses Lennard-Jones units \n" << std::endl;
  std::cout << "Reading parameters from " << simParameters << std::endl;

  // Read seed & primes for random numbers
  int p1, p2;
  std::ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  std::ifstream Seed("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];

  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;
  std::cout << "Temperature = " << temp << std::endl;

  ReadInput >> npart;
  std::cout << "Number of particles = " << npart << std::endl;

  ReadInput >> rho;
  std::cout << "Density of particles = " << rho << std::endl;
  vol = (double)npart/rho;
  box = std::pow(vol,1.0/3.0);
  std::cout << "Volume of the simulation box = " << vol << std::endl;
  std::cout << "Edge of the simulation box = " << box << std::endl;

  ReadInput >> rcut;
  std::cout << "Cutoff of the interatomic potential = " << rcut << std::endl << std::endl;

  //Tail corrections for potential energy and pressure
  vtail = (8.0*M_PI*rho)/(9.0*std::pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*std::pow(rcut,3));
  ptail = (32.0*M_PI*rho)/(9.0*std::pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*std::pow(rcut,3));
  std::cout << "Tail correction for the potential energy = " << vtail << std::endl;
  std::cout << "Tail correction for the pressure           = " << ptail << std::endl;

  ReadInput >> delta;

  ReadInput >> nstep;

  ReadInput >> nblk;

  std::cout << "The program perform Metropolis moves with uniform translations" << std::endl;
  std::cout << "Moves parameter = " << delta << std::endl;
  std::cout << "Number of blocks = " << nblk << std::endl;
  std::cout << "Number of steps in one block = " << nstep << std::endl << std::endl;
  ReadInput.close();

  for (auto& el : props) {
    walker[el] = 0.;
    block_average[el] = 0.;
    glob_average[el] = 0.;
    glob_average2[el] = 0.;
  }

  //Read initial configuration
  std::cout << "Read initial configuration from file "
  << initial_configuration << std::endl << std::endl;

  double tmpx, tmpy, tmpz;
  x.resize(npart);
  y.resize(npart);
  z.resize(npart);
  for (unsigned int i=0; i<npart; ++i){
    ReadConf >> tmpx >> tmpy >> tmpz;
    x[i] = Pbc( tmpx * box );
    y[i] = Pbc( tmpy * box );
    z[i] = Pbc( tmpz * box );
  }
  ReadConf.close();

  bin_size = (box/2.0)/(double)nbins;

  std::ofstream binning("results/binning.dat");
  for(unsigned int i=0; i<nbins; ++i) {
    binning << i*bin_size+bin_size/2. << std::endl;
  }
  binning.close();

  histo_walker.resize(nbins);
  histo_block_average.resize(nbins);
  histo_glob_average.resize(nbins);
  histo_glob_average2.resize(nbins);

  //Evaluate potential energy and pressure of the initial configuration
  Measure();
  //Print initial values for the potential energy and pressure
  std::cout << "Initial potential energy (with tail corrections) = "
  << walker.at("energy")/(double)npart + vtail << std::endl;
  std::cout << "Virial                   (with tail corrections) = "
  << walker.at("pressure")/(double)npart + ptail << std::endl;
  std::cout << "Pressure                 (with tail corrections) = "
  << rho*temp+(walker.at("pressure")+(double)npart*ptail)/vol
  << std::endl << std::endl;


  blk_norm=0;
  attempted=0;
  accepted=0;
  Epot.open("results/epot.dat");
  Pres.open("results/pres.dat");
  Gerr.open("results/gerr.dat");
  Gave.open("results/gave.dat");

  if(instant){
    ist_pot.open("results/instant_epot.dat");
    ist_pres.open("results/instant_pres.dat");
  }
}

double MolDynMC::Boltzmann(double xx, double yy, double zz, unsigned int ip) {
  double ene=0.0;
  double dx, dy, dz, dr;

  for(unsigned int  i=0; i<npart; ++i)
  {
    if(i != ip)
    {
      // distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = std::sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/std::pow(dr,12) - 1.0/std::pow(dr,6);
      }
    }
  }
  return 4.0*ene;
}

void MolDynMC::Move(void) {
  int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;


  for(unsigned int  i=0; i<npart; ++i) {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (int)(rnd.Rannyu()*npart);

    //Old
    xold = x[o];
    yold = y[o];
    zold = z[o];

    energy_old = Boltzmann(xold,yold,zold,o);

    //New
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

    energy_new = Boltzmann(xnew,ynew,znew,o);

    //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu())
    {
      //Update
      x[o] = xnew;
      y[o] = ynew;
      z[o] = znew;

      accepted += 1.0;
    }
    attempted += 1.0;
  }
}

void MolDynMC::ConfFinal() {
  std::ofstream WriteConf;

  std::cout << "Print final configuration to file config.final " << std::endl << std::endl;
  WriteConf.open("config.final");
  for(unsigned int  i=0; i<npart; ++i) {
    WriteConf << x[i]/box << " " <<  y[i]/box << " " << z[i]/box << std::endl;
  }
  WriteConf.close();
  rnd.SaveSeed();
}

void MolDynMC::ConfXYZ(unsigned int nconf) { //Write configuration in .xyz format
  std::ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + std::to_string(nconf) + ".xyz");
  WriteXYZ << npart << std::endl;
  WriteXYZ << "This is only a comment!" << std::endl;
  for(unsigned int  i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << std::endl;
  }
  WriteXYZ.close();
}

void MolDynMC::Accumulate(void) {
  //Update block averages
  for(auto& el: props) {
    block_average.at(el) = block_average.at(el) + walker.at(el);
  }
  for(unsigned int i=0; i<nbins; ++i) {
    histo_block_average[i] += histo_walker[i];
  }
  blk_norm ++;
}

void MolDynMC::Reset(unsigned int iblk) {
  //Reset block averages
  // we choose to lookup the maps' elements with .at for safety:
  // If you access a key using the indexing operator [] that is
  // not currently a part of a map, then it automatically adds
  // a key for you!!
  if(iblk == 0) {
    for(auto & el : props) {
      glob_average.at(el) = 0.;
      glob_average2.at(el) = 0.;
    }
    for(auto& el : histo_glob_average) {
      el=0.;
    }
    for(auto& el : histo_glob_average2) {
      el=0.;
    }
  }

  for(auto & el : props) block_average.at(el) = 0.;

  for(auto & el : histo_block_average) el = 0.;

  blk_norm = 0.;
  attempted = 0.;
  accepted = 0.;
}

void MolDynMC::Run() {
  unsigned int nconf = 1;
  //Simulation
  for(unsigned int iblk=1; iblk <= nblk; ++iblk) {
    Reset(iblk);   //Reset block averages
    for(unsigned int istep=1; istep <= nstep; ++istep) {
      Move();
      Measure(istep);
      Accumulate(); //Update block averages

      if(istep%10 == 0) {
        ConfXYZ(nconf);
        nconf++;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
}

void MolDynMC::Averages(unsigned int iblk) {
  //Print results for current block
  // double r, gdir;

  std::cout << "Block number " << iblk << std::endl;
  std::cout << "Acceptance rate " << (double)accepted/attempted << std::endl
            << std::endl;

  //Potential energy
  double stima_pot, err_pot;
  stima_pot = block_average.at("energy")/blk_norm/(double)npart + vtail;
  glob_average.at("energy") += stima_pot;
  glob_average2.at("energy") += stima_pot*stima_pot;
  err_pot = Error(glob_average.at("energy"), glob_average2.at("energy"), iblk);

  //Pressure
  double stima_pres, err_pres;
  stima_pres = rho * temp + (block_average.at("pressure")/blk_norm + ptail * (double)npart) / vol;
  glob_average.at("pressure") += stima_pres;
  glob_average2.at("pressure") += stima_pres*stima_pres;
  err_pres = Error(glob_average.at("pressure"), glob_average.at("pressure"), iblk);

  //Potential energy per particle
  Epot << " " << iblk <<  " " << stima_pot << " "
       << glob_average.at("pressure")/(double)iblk << " " << err_pot
       << std::endl;

  //Pressure
  Pres << " " << iblk <<  " " << stima_pres << " "
  << glob_average.at("pressure")/(double)iblk << " " << err_pres
  << std::endl;

  // err_g
  double stima_g; //, err_g;
    std::vector<double> err_g;
    for(unsigned int i=0; i<nbins; ++i){
        stima_g = (double)histo_block_average[i]/blk_norm;
        histo_glob_average[i] += stima_g;
        histo_glob_average2[i] += stima_g*stima_g;
        err_g.push_back(Error(histo_glob_average[i], histo_glob_average2[i],iblk));
    }
    Gerr << std::endl;

  //g(r)
  double norm = 0.; // *(4.*M_PI/3.);
  for(auto & it : histo_glob_average) norm += it;

  for(auto & it : histo_glob_average) Gave << it/(norm) << " ";
  Gave << std::endl;

  for(auto & it : err_g) Gerr << it/(norm)<< " ";
  Gerr << std::endl;

  std::cout << "----------------------------" << std::endl << std::endl;

}

void MolDynMC::Measure(unsigned int istep) {
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

  //reset the hystogram of g(r)
  for(auto &el: histo_walker) el=0;

  //cycle over pairs of particles
  for (unsigned int i=0; i<npart-1; ++i) {
    for (unsigned int j=i+1; j<npart; ++j) {
      // distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      //update of the histogram of g(r)
         for (unsigned int k=0; k<nbins; ++k){
             if(dr>bin_size*k and dr<=bin_size*(k+1)){
//                  r=std::sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
                 histo_walker[k]+=3./(std::pow(dr+bin_size,3)-std::pow(dr,3))
                     /(2.*M_PI);
                 break;
             }
         }


      if(dr < rcut) {
        vij = 1.0/std::pow(dr,12) - 1.0/std::pow(dr,6);
        wij = 1.0/std::pow(dr,12) - 0.5/std::pow(dr,6);
        // contribution to energy and pressure
        v += vij;
        w += wij;
      }
    }
  }
  walker.at("energy") = 4.0 * v;
  walker.at("pressure") = 48.0 * w / 3.0;

  if(instant) {
    double ener, pre;
    ener = walker.at("energy")/(double)npart + vtail;
    pre = rho*temp+(walker.at("pressure")+ ptail*(double)npart)/vol;
    ist_pot << istep << " " << ener << std::endl;
    ist_pres << istep << " " << pre << std::endl;
  }
}

MolDynMC::~MolDynMC(){
    Epot.close();
    Pres.close();
    Gerr.close();
    Gave.close();
    if(instant) {
        ist_pot.close();
        ist_pres.close();
    }
}
