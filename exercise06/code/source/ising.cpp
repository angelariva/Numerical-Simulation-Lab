#include "ising.h"

Ising1D::Ising1D(std::string old_configuration) {
  rnd.SetRandom("Primes", "seed.in");
  props = {"energy", "capacity", "magnetization", "susceptibility"};
  for (auto& el : props) {
    walker[el] = 0.;
    block_average[el] = 0.;
    glob_average[el] = 0.;
    glob_average2[el] = 0.;
  }

  Input();

  if(metro) {
    Move = &Ising1D::MetropolisMove;
    std::cout << "Metropolis sampling " << std::endl;
  } else {
    Move = &Ising1D::GibbsMove;
    std::cout << "Gibbs sampling " << std::endl;
  }

  if(old_configuration=="") { //start from random config
/*    Ene.open("outputs/ene.dat");   //delete content of file
    Heat.open("outputs/heat.dat");
    Mag.open("outputs/mag.dat");
    Chi.open("outputs/chi.dat");*/
    //initial configuration
    s.resize(n_spin);
    for(auto& it : s){
        if(rnd.Rannyu() >= 0.5) it=1;
        else it = -1;
    }
  }else{
/*    Ene.open("outputs/ene.dat", std::ios::app);
    Heat.open("outputs/heat.dat", std::ios::app);
    Mag.open("outputs/mag.dat", std::ios::app);
    Chi.open("outputs/chi.dat", std::ios::app);*/
    std::ifstream input_file(old_configuration);
    if(input_file.fail()){
          std::cout << "Unable to open " << old_configuration << "\n";
          exit(1);
    }
    double tmp;
    input_file >> tmp;
    while(!input_file.eof()){
        s.push_back(tmp);
        input_file >> tmp;
    }
    input_file.close();
    if (s.size()!=n_spin){
        std::cout << "Spin configurations found  in '"
                  << old_configuration
                  << "' don't belong to a simulation compatible with "
                  << "parameters defined in 'input.dat'\n";
        exit(0);
    }
  }

  //Evaluate energy etc. of the initial configuration
  Measure();

  //Print initial values for the potential energy and virial
  std::cout << "Initial energy = " << walker["energy"]/(double)n_spin << std::endl;

}

Ising1D::~Ising1D() {
  Ene.close();
  Heat.close();
  Mag.close();
  Chi.close();
}


void Ising1D::Input()  {
  std::ifstream ReadInput;
  if(ReadInput.fail()){
    std::cout << "Unable to open 'input.dat'" << "\n";
    exit(1);
  }

  std::cout << "Classic 1D Ising model             " << std::endl;
  std::cout << "Monte Carlo simulation             " << std::endl << std::endl;
  std::cout << "Nearest neighbour interaction      " << std::endl << std::endl;
  std::cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << std::endl << std::endl;
  std::cout << "The program uses k_B=1 and mu_B=1 units " << std::endl;

//Read seed for random numbers
   int p1, p2;
   std::ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   std::ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  std::cout << "Temperature = " << temp << std::endl;

  ReadInput >> n_spin;
  std::cout << "Number of spins = " << n_spin << std::endl;

  ReadInput >> J;
  std::cout << "Exchange interaction = " << J << std::endl;

  ReadInput >> h;
  std::cout << "External field = " << h << std::endl << std::endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs
  ReadInput >> nblk;
  ReadInput >> nstep;

  std::cout << "Number of blocks = " << nblk << std::endl;
  std::cout << "Number of steps in one block = " << nstep << std::endl << std::endl;
  ReadInput.close();

}


void Ising1D::Measure(int step)
{
  // int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (unsigned int i=0; i<n_spin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker.at("energy") = u;
  walker.at("capacity") = u*u;
  walker.at("magnetization") = m;
  walker.at("susceptibility") = m*m;

  if(step !=0) {
    std::ofstream EneIn, CapIn, MagIn, SusIn;
    EneIn.open("results/instant.ene.0",std::ios::app);
    CapIn.open("results/instant.cap.0",std::ios::app);
    MagIn.open("results/instant.mag.0",std::ios::app);
    SusIn.open("results/instant.sus.0",std::ios::app);

    EneIn << step << " " << walker.at("energy") << std::endl;
    CapIn << step << " " << walker.at("capacity") << std::endl;
    MagIn << step << " " << walker.at("magnetization") << std::endl;
    SusIn << step << " " << walker.at("susceptibility") << std::endl;

    EneIn.close();
    CapIn.close();
    MagIn.close();
    SusIn.close();
  }
}

void Ising1D::Reset(unsigned int iblk) {
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
   }

   for(auto & el : props) block_average.at(el) = 0.;

   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Ising1D::Accumulate(void) {
//Update block averages

   for(auto& el: props)
    block_average.at(el) = block_average.at(el) + walker.at(el);

   blk_norm ++;
}

void Ising1D::BlockAverages(unsigned int iblk) {
  //Print results for current block
  std::ofstream Ene, Heat, Mag, Chi;

  std::cout << "Block number " << iblk << std::endl;
  std::cout << "Acceptance rate " << accepted/attempted << std::endl << std::endl;

  Ene.open("results/output.ene.0",std::ios::app);  // Energy
  stima_u = block_average.at("energy")/blk_norm/(double)n_spin;
  glob_average.at("energy")  += stima_u;
  glob_average2.at("energy") += stima_u*stima_u;
  err_u = Error(glob_average.at("energy"), glob_average2.at("energy"), iblk);
  Ene << " " << iblk <<  " " << stima_u << " "
  << glob_average.at("energy")/(double)iblk << " " << err_u << std::endl;
  Ene.close();

  Heat.open("results/output.heat.0",std::ios::app);  // Heat capacity
  stima_c = beta*beta*(block_average.at("capacity")/blk_norm);
  stima_c = (stima_c-stima_u*stima_u*(double)n_spin)/(temp*temp);
  glob_average.at("capacity")  += stima_c;
  glob_average2.at("capacity") += stima_c*stima_c;
  err_c = Error(glob_average.at("capacity"), glob_average2.at("capacity"), iblk);
  Heat << " " << iblk <<  " " << stima_c << " "
  << glob_average.at("capacity")/(double)iblk << " " << err_c << std::endl;
  Heat.close();

  Mag.open("results/output.mag.0",std::ios::app);  // Magnetization
  stima_m = block_average["magnetization"]/blk_norm/(double)n_spin;
  glob_average.at("magnetization")  += stima_m;
  glob_average2.at("magnetization") += stima_m*stima_m;
  err_m = Error(glob_average.at("magnetization"), glob_average2.at("magnetization"), iblk);
  Mag << " " << iblk <<  " " << stima_m << " " << glob_average.at("magnetization")/(double)iblk
  << " " << err_m << std::endl;
  Mag.close();

  Chi.open("results/output.chi.0",std::ios::app);  // susceptibility
  stima_x = beta*block_average.at("susceptibility")/blk_norm/(double)n_spin;
  glob_average.at("susceptibility")  += stima_x;
  glob_average2.at("susceptibility") += stima_x*stima_x;
  err_x = Error(glob_average.at("susceptibility"), glob_average2.at("susceptibility"), iblk);
  Chi << " " << iblk <<  " " << stima_x << " " << glob_average.at("susceptibility")/(double)iblk
  << " " << err_x << std::endl;
  Chi.close();
}

void Ising1D::ConfFinal(void) {
  std::ofstream WriteConf;

  std::cout << "Print final configuration to file config.final " << std::endl << std::endl;
  WriteConf.open("config.final");
  for (auto& el : s)
    WriteConf << el << std::endl;

  WriteConf.close();
  rnd.SaveSeed();
}

void Ising1D::MetropolisMove() {
  int flip;
  double deltaen, Q;

  //Select randomly a particle (for C++ syntax, 0 <= flip <= n_spin-1)
  flip = (int)(rnd.Rannyu()*n_spin);
  deltaen = -2*Boltzmann(s[flip], flip);
  Q = std::exp(-beta*deltaen);
  double A = std::min(1., Q);
  if(rnd.Rannyu() <= A) {
    s[flip] = -s[flip];
    accepted++;
  }
  attempted++;
}

void Ising1D::GibbsMove() {
  int flip;
  double p, deltaen1, deltaen2, Q;
  //Select randomly a particle (for C++ syntax, 0 <= o <= n_spin-1)
  flip = (int)(rnd.Rannyu()*n_spin);
  deltaen1 = -2.*Boltzmann(1, flip);
  deltaen2 = -2.*Boltzmann(-1,flip);
  Q = std::exp(-beta*deltaen1)  + std::exp(-beta*deltaen2);
  p = 1./(1.+Q);
  if(rnd.Rannyu() <= p){
    s[flip] = 1.;
    accepted++;
  } else s[flip] = -1.;
  attempted++;
}

void Ising1D::Run(bool instant) {
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      (*this.*Move)();
      if(instant) Measure(istep);
      else Measure();
      Accumulate(); //Update block averages
    }
    BlockAverages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
}
