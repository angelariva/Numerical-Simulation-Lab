#include "MolDyn_NVE.h"

using namespace std;

MolecularDynamics::MolecularDynamics(std::string simParameters, std::string configFile) {

  Input(simParameters);
  // Read initial configuration from the configuration file
  // we are starting from scratch our simulation in this CONSTRUCTOR
  cout << "Read initial configuration from file " + configFile << endl << endl;

  std::ifstream ReadConf(configFile);
  if(ReadConf.fail()){
      std::cerr << "Error opening config file\n";
      exit(1);
  }
  for (unsigned int i=0; i<npart; ++i){
       ReadConf >> x[i] >> y[i] >> z[i];
       x[i] = x[i] * box;
       y[i] = y[i] * box;
       z[i] = z[i] * box;
   }
   ReadConf.close();

   //Prepare initial velocities
   cout << "Prepare random velocities with Maxwell-Boltzmann distribution" << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (unsigned int i=0; i<npart; ++i){
       vx[i] = rand.Rannyu() - 0.5;
       vy[i] = rand.Rannyu() - 0.5;
       vz[i] = rand.Rannyu() - 0.5;

       sumv[0] += vx[i];
       sumv[1] += vy[i];
       sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (unsigned int i=0; i<npart; ++i){
       vx[i] = vx[i] - sumv[0];
       vy[i] = vy[i] - sumv[1];
       vz[i] = vz[i] - sumv[2];

       sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * double(temp) / sumv2);   // fs = velocity scale factor

   for (unsigned int i=0; i<npart; ++i) {
       vx[i] *= fs;
       vy[i] *= fs;
       vz[i] *= fs;
       xold[i] = Pbc(x[i] - vx[i] * delta);
       yold[i] = Pbc(y[i] - vy[i] * delta);
       zold[i] = Pbc(z[i] - vz[i] * delta);
   }

}

MolecularDynamics::MolecularDynamics(std::string simParameters, std::string configFile, std::string oldConfigFile) {
  Input(simParameters);
  // Read initial configuration from the configuration file (final in this case!)
  cout << "Read initial configuration from file " + configFile << endl << endl;

  std::ifstream ReadConf(configFile);
  if(ReadConf.fail()){
      std::cerr << "Error opening config file\n";
      exit(1);
  }
  for (unsigned int i=0; i<npart; ++i){
       ReadConf >> x[i] >> y[i] >> z[i];
       x[i] = x[i] * box;
       y[i] = y[i] * box;
       z[i] = z[i] * box;
   }
   ReadConf.close();
   ReadConf.clear();

   // read previous configuration file (from old simulation)
   cout << "Read old configuration from file " + oldConfigFile << endl << endl;
   ReadConf.open(oldConfigFile);
   if(ReadConf.fail()){
      std::cerr << "Error opening config file\n";
      exit(1);
   }
   for (unsigned int i=0; i<npart; ++i){
         ReadConf >> xold[i] >> yold[i] >> zold[i];
         xold[i] = xold[i] * box;
         yold[i] = yold[i] * box;
         zold[i] = zold[i] * box;
   }
   ReadConf.close();
   ReadConf.clear();

   Move();   // we compute the forces acting on the particles to use them in
              // the Verlet algorithm (take a look at the loop below:)
              // x(t+dt) = 2x(t) - x(t-dt) + F(t)dt^2/m (with PBCs)

   // ESTIMATE for TEMPERATURE
   double t = 0.;
   for (unsigned int i=0; i<npart; ++i)
       t += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   stima_temp = t/(npart*3.);

   double fs = sqrt(temp / stima_temp);   // fs = velocity scale factor
   cout << "scaling factor = " << fs << endl;
   // upadating velocities and positions
   for (unsigned int i=0; i<npart; ++i) {
       vx[i] *= fs;
       vy[i] *= fs;
       vz[i] *= fs;

       double xi = xold[i];
       double yi = yold[i];
       double zi = zold[i];

       xold[i] = Pbc(x[i] - 2.*vx[i] * delta);
       yold[i] = Pbc(y[i] - 2.*vy[i] * delta);
       zold[i] = Pbc(z[i] - 2.*vz[i] * delta);

       x[i] = xi;
       y[i] = yi;
       z[i] = zi;

   }

}

MolecularDynamics::~MolecularDynamics() {
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Press.close();
}

void MolecularDynamics::Input(std::string simParameters) {
  Epot.open("results/output_epot.dat");
  Ekin.open("results/output_ekin.dat");
  Temp.open("results/output_temp.dat");
  Etot.open("results/output_etot.dat");
  Press.open("results/output_press.dat");

  rand.SetRandom("Primes", "seed.in");
  std::ifstream ReadInput(simParameters);
  if(ReadInput.fail()){
      std::cerr << "Unable to open setup file\n";
      exit(1);
  }

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  // reading simulation & physical parameters from the input file
  ReadInput >> temp;
  ReadInput >> npart;
  // resizing of vectors containing simulations' info
  x.resize(npart);
  y.resize(npart);
  z.resize(npart);
  xold.resize(npart);
  yold.resize(npart);
  zold.resize(npart);
  vx.resize(npart);
  vy.resize(npart);
  vz.resize(npart);
  fx.resize(npart);
  fy.resize(npart);
  fz.resize(npart);

  cout << "Number of particles = " << npart << endl;
  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;

  ReadInput >> measure_time_interval;
  cout << "Measures performed every " << measure_time_interval <<
  " time steps." << endl;
  unsigned int n_blocks;
  ReadInput >> n_blocks;
  std::cout <<"Statistical error computed with "<< n_blocks << " blocks." << endl << endl;

  // resizing of vectors to perform the blocking averages
  est_pot.resize(n_blocks);
  est_kin.resize(n_blocks);
  est_etot.resize(n_blocks);
  est_temp.resize(n_blocks);
  est_press.resize(n_blocks);
  // initializing them to zero
  std::fill(est_pot.begin(), est_pot.end(), 0.);
  std::fill(est_kin.begin(), est_kin.end(), 0.);
  std::fill(est_etot.begin(), est_etot.end(), 0.);
  std::fill(est_temp.begin(), est_temp.end(), 0.);
  std::fill(est_press.begin(), est_press.end(), 0.);

  ReadInput.close();

  // compute block size (for the blocking averages) and setting block index
  // (L = M/N of the previous exercises, M is nstep/measure_time_interval, N is
  // nblocks)
  block_size=(nstep/measure_time_interval)/n_blocks;
  iblock=0;
  imeasure=0;
  bin_size = (box/2.0)/(double)nbins;

  Binn.open("results/binning.dat");

  histo.resize(nbins);
  for(auto & el : histo) {
    el.resize(n_blocks);
    std::fill(el.begin(), el.end(), 0.);
  }
}


void MolecularDynamics::ConfXYZ(int nconf) const { //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (unsigned int i=0; i<npart; ++i)
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  WriteXYZ.close();
}

void MolecularDynamics::ConfFinal(std::string filename) const {
  ofstream WriteConf;
  cout << "Print configuration to " + filename << endl;
  WriteConf.open(filename);

  for(unsigned int i=0; i<npart; i++) WriteConf << x[i]/box << " " << y[i]/box << " " << z[i]/box << endl;
  WriteConf.close();
}

void MolecularDynamics::Move(){ //Move particles with Verlet algorithm
  double xnew, ynew, znew;
  Force();

  for(unsigned int i=0; i<npart; ++i){ //Verlet integration scheme
    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
}

void MolecularDynamics::Force() {
  // Compute forces as -Grad_ip V(r)
  double d1, d2, d3, dr;
  double multiplier;
  fill(fx.begin(), fx.end(), 0.);
  fill(fy.begin(), fy.end(), 0.);
  fill(fz.begin(), fz.end(), 0.);

  for(unsigned int j=0; j<npart; ++j) {
    for(unsigned int i=0; i<npart; ++i) {
      if(i != j) {
        d1 = Pbc( x[j] - x[i] );  // distance j-i in pbc
        d2 = Pbc( y[j] - y[i] );
        d3 = Pbc( z[j] - z[i] );

        dr = d1*d1 + d2*d2 + d3*d3;
        dr = sqrt(dr);

        if(dr < rcut){
          // -Grad_ip V(r)
          multiplier = (48.0/pow(dr,14) - 24.0/pow(dr,8));
          fx[j] += d1*multiplier;
          fy[j] += d2*multiplier;
          fz[j] += d3*multiplier;
        }
      }
    }
  }
}

void MolecularDynamics::Measure(){
  //Properties measurement
  double v, t;
  double dx, dy, dz, dr, r;

  v = 0.0; //reset observables
  t = 0.0;
  stima_press = 0.;
  //cycle over pairs of particles
  for (unsigned int i=0; i<npart-1; ++i) {
    for (unsigned int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     //update of the histogram of g(r)
     for (unsigned int k=0; k<nbins; ++k){
       if(dr>bin_size*k and dr<bin_size*(k+1)){
         r=std::sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
         histo[k][iblock]+=1./(std::pow(r+dr,3)-std::pow(r,3));
         break;
       }
     }

     if(dr < rcut){
       // Potential energy
       v += 4.0/pow(dr,12) - 4.0/pow(dr,6);
       stima_press += 1./pow(dr,12)-0.5/pow(dr,6);
     }
    }
  }

  //Kinetic energy
  for(unsigned int i=0; i<npart; ++i)
    t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
  stima_press = 16.*stima_press/(vol);
  stima_press += stima_temp*rho;

  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  Press << stima_press << endl;

  imeasure++;
  iblock = imeasure/block_size;
  est_pot[iblock] += stima_pot;
  est_kin[iblock] += stima_kin;
  est_temp[iblock] += stima_temp;
  est_etot[iblock] += stima_etot;
  est_press[iblock] += stima_press;
}

void MolecularDynamics::BlockingResults() {
  for(auto &el : est_pot) el /= block_size;
  for(auto &el : est_kin) el /= block_size;
  for(auto &el : est_temp) el /= block_size;
  for(auto &el : est_etot) el /= block_size;
  for(auto &el : est_press) el /= block_size;
  for(auto &el : histo) {
    for(auto & k : el) {
      k /= block_size;
    }
  }

  vector<double> pot_err = blocking_error(est_pot);
  vector<double> kin_err = blocking_error(est_kin);
  vector<double> temp_err = blocking_error(est_temp);
  vector<double> etot_err = blocking_error(est_etot);
  vector<double> press_err = blocking_error(est_press);
  std::vector<std::vector<double> > histo_err(histo.size());
  for (unsigned int i=0; i<histo.size(); ++i) {
    histo_err[i]=blocking_error(histo[i]);
  }

  Gerr.open("results/gerr.dat"),
  Gave.open("results/gave.dat");

  for(unsigned int i = 0; i<histo[0].size(); ++i) {
    for(unsigned int j = 0; j<histo.size(); ++j) {
      Gave << histo[j][i] << " ";
      Gerr << histo_err[j][i] << " ";
    }
    Gave << endl;
    Gerr << endl;
  }

  for(unsigned int i=0; i<nbins; ++i) {
    double a;
    a = i*bin_size+bin_size/2.;
    Binn << a << " ";
  }
  Binn << endl;
  Gerr.close();
  Gave.close();

  ofstream out("results/ave_epot.dat");
  for(unsigned int i=0; i<est_pot.size(); ++i) {
    out << i << " " << est_pot[i] << " " << pot_err[i] << endl;
  }
  out.close();

  out.open("results/ave_ekin.dat");
  for(unsigned int i=0; i<est_kin.size(); ++i) {
     out << i << " " << est_kin[i] << " " << kin_err[i] << endl;
  }
  out.close();

  out.open("results/ave_etot.dat");
  for(unsigned int i=0; i<est_etot.size(); ++i) {
     out << i << " " << est_etot[i] << " " << etot_err[i] << endl;
  }
  out.close();

  out.open("results/ave_temp.dat");
  for(unsigned int i=0; i<est_temp.size(); ++i) {
     out << i << " " << est_temp[i] << " " << temp_err[i] << endl;
  }
  out.close();

  out.open("results/ave_press.dat");
  for(unsigned int i=0; i<est_press.size(); ++i) {
     out << i << " " << est_press[i] << " " << press_err[i] << endl;
   }
  out.close();

}

void MolecularDynamics::RunSimulation() {
  unsigned int nconf=1; // starting from configuration number 1
  std::ofstream out("frames/config_1.xyz"); //check if frames folder exists
  if(out.fail()){
      std::cout << "\n\nERROR: Unable to write on folder 'frames', create it"
          "and run again!\n\n";
        exit(2);
  }
  out.close();
  // modified in order to print one xyz file less: the corresponding
  // configuration will be printed in the last lines to "old.0" and
  // "config.final".
  for (unsigned int istep=1; istep<nstep; ++istep) {
      Move();
      if(istep % iprint ==0){
          cout << "Number of time-steps: " << istep << endl;
          ConfXYZ(nconf);
        }
      if(istep % measure_time_interval == 0) {
            Measure();
            nconf++;
        }
  }
  cout << "Number of time-steps: " << nstep << endl << endl;
  ConfFinal("old.0");
  Move();
  Measure();
  ConfFinal("config.final");
  BlockingResults();
}

double MolecularDynamics::Pbc(double r) const {
    return r - box * rint(r/box);
}
