/*
Angela Riva
EXERCISE 05
We use the metropolis algorithm to sample the probabilites
for the electron in the hidrogen atom of finding itself in:
- 1s ground state
- 2p state
To do that we use as a transition probability T(x|y) a gaussian
and a uniform distribution. We calculate <r> in the two states
(with data blocking).
The scripts will help us in studying how the starting point
choice reflects in the equilibration phase of the simulation.
The quantities are rescaled in Bohr radius units
*/

#include "hydrogen.h"

void help_message(char* argv[]){
    std::cout <<"USAGE: \t" << argv[0] << " [option1] [options2] \n\n";
    std::cout << "EXAMPLE: ./main --pos 1.22 3 0.7 2 \n\n";
    std::cout << "OPTION 1:" << std::endl;
    std::cout <<"\t--help\t\t\tto print out this help message\n" <<std::endl;
    std::cout <<"\t--pos \t\t\tto print on file the positions sampled"<<std::endl;
    std::cout <<"\t--rad \t\t\tto print on file the average radius"<<std::endl;
    std::cout <<"\t--pos+rad \t\tto print on file the pos sampled & av radius"<<std::endl;
    std::cout <<"\t--far \t\t\tto print on file the pos sampled from a position\n";
    std::cout <<"\t\t\t\tfar, far away from the origin\n"<<std::endl;
    std::cout << "OPTIONs 2, in order: [step100unif] [step210unif] [step100gaus] [step210gaus]" << std::endl;
    std::cout <<"\t Double values to set the step of: "<<std::endl;
    std::cout <<"\t [step100unif]: ground state distribution with uniform transition probability; "<<std::endl;
    std::cout <<"\t [step210unif]: first excited state distribution with uniform transition probability;"<<std::endl;
    std::cout <<"\t [step100gaus]: ground state distribution with gaussian transition probability; "<<std::endl;
    std::cout <<"\t [step210gaus]: first excited state distribution with gaussian transition probability."<<std::endl;
    exit(0);
}


int main(int argc, char *argv[]) {
  //Options setup
  if(argc!=6){
      if(argc==2 and std::string(argv[1])=="--help") help_message(argv);
      else {
        std::cout << "Unrecognized option \n" << std::endl;
        help_message(argv);
      }
  }
  else{
    if (std::string(argv[1])=="--help")
      help_message(argv);
    if ((std::string(argv[1]))!="--pos" and (std::string(argv[1]))!="--rad" and
        (std::string(argv[1]))!="--pos+rad" and (std::string(argv[1]))!="--far") {
      std::cout << "Unrecognized option" << std::endl;
      help_message(argv);
    }
  }

  unsigned int M = 1e7;
  unsigned int N = 1e2;
  unsigned int L = int(M/N);

  double step100uni = atof(argv[2]);
  double step210uni = atof(argv[3]);
  double step100gaus = atof(argv[4]);
  double step210gaus = atof(argv[5]);

  std::cout << "step 100 uniform = " << step100uni << std::endl;
  std::cout << "step 210 uniform = " << step210uni << std::endl;
  std::cout << "step 100 gaussian = " << step100gaus << std::endl;
  std::cout << "step 210 gaussian = " << step210gaus << std::endl;

  double x=0.,y=0.,z=0.;


  if(std::string(argv[1])=="--far") {
    x = 100., y=100., z=100.;
    M = 10000;

    std::cout << "\nPrinting the sampled radius, starting far, far away (100, 100, 100)\n";
    std::cout << "from the origin and using UNIFORM transition probability: " <<std::endl;
    hydrogen hydro(x,y,z,step100uni, "100", "uniform");
    std::vector<double> r100(M, 0.), r210(M, 0.);

    for (unsigned int i=0; i< M; ++i){
      hydro.metropolis();
      r100[i] = hydro.get_radius();
    }
    std::cout << "Uniform 100 acceptance: " << hydro.acceptance() << std::endl;

    hydro.reset_metropolis(x,y,z,step210uni, "210", "uniform");
    for (unsigned int i=0; i<M; ++i){
      hydro.metropolis();
      r210[i] = hydro.get_radius();
    }
    std::cout << "Uniform 210 acceptance: " << hydro.acceptance() << std::endl;

    std::ofstream out("radius100far.dat");
    for(unsigned int i=0; i<M; ++i)
      out << i << " " << r100[i] << std::endl;
    out.close();

    out.open("radius210far.dat");
    for(unsigned int i=0; i<M; ++i)
      out << i << " " << r210[i] << std::endl;
    out.close();
  }

  if(std::string(argv[1])=="--rad" or std::string(argv[1])=="--pos+rad") {
    std::cout << "\nPrinting the average radius, starting from the origin (0,0,0) \n";
    std::cout << "and using UNIFORM transition probability: "<<std::endl;

    hydrogen hydro(x,y,z,step100uni, "100", "uniform");
    std::vector<double> r100(N, 0.), r210(N, 0.);

    for (unsigned int i=0; i< N; ++i){
      for(unsigned j=0; j<L; ++j){
        hydro.metropolis();
        r100[i] += hydro.get_radius();
      }
      r100[i]/=L;
    }
    std::cout << "Uniform 100 acceptance: " << hydro.acceptance() << std::endl;

    hydro.reset_metropolis(x,y,z,step210uni, "210", "uniform");
    for (unsigned int i=0; i< N; ++i){
      for(unsigned j=0; j<L; ++j){
         hydro.metropolis();
         r210[i]+=hydro.get_radius();
       }
       r210[i]/=L;
    }
    std::cout << "Uniform 210 acceptance: " << hydro.acceptance() << std::endl;
    //computing errors with blocking method
    std::vector<double> error100=blocking_error(r100);
    std::vector<double> error210=blocking_error(r210);

    std::ofstream out("radius100uniform.dat");
    for(unsigned int i=0; i<N; ++i)
      out << r100[i] << " " << error100[i] << std::endl;
    out.close();

    out.open("radius210uniform.dat");
    for(unsigned int i=0; i<N; ++i)
      out << r210[i] << " " << error210[i] << std::endl;
    out.close();

    /* Mean radius computation with gaussian random numbers */
    std::fill(r100.begin(), r100.end(), 0.);
    std::fill(r210.begin(), r210.end(), 0.);

    hydro.reset_metropolis(x,y,z,step100gaus,"100", "gaussian");
    std::cout << "\nPrinting the average radius, starting from the origin (0,0,0)\n";
    std::cout << "and using GAUSSIAN transition probability: "<<std::endl;

    for (unsigned int i=0; i< N; ++i) {
      for(unsigned j=0; j<L; ++j) {
        hydro.metropolis();
        r100[i]+=hydro.get_radius();
      }
      r100[i]/=L;
    }
    std::cout << "Gauss 100 acceptance: " << hydro.acceptance() << std::endl;

    hydro.reset_metropolis(x,y,z,step210gaus, "210", "gaussian");
    for (unsigned int i=0; i< N; ++i) {
      for(unsigned j=0; j<L; ++j) {
        hydro.metropolis();
        r210[i]+=hydro.get_radius();
      }
      r210[i]/=L;
    }
    std::cout << "Gauss 210 acceptance: " << hydro.acceptance() << std::endl;
    error100=blocking_error(r100);
    error210=blocking_error(r210);

    out.open("radius100gauss.dat");
    for(unsigned int i=0; i<N; ++i)
      out << r100[i] << " " << error100[i] << std::endl;
    out.close();

    out.open("radius210gauss.dat");
    for(unsigned int i=0; i<N; ++i)
      out << r210[i] << " " << error210[i] << std::endl;
    out.close();
  }

  if(std::string(argv[1])=="--pos" or std::string(argv[1])=="--pos+rad") {

    N=1e5;
    std::cout << "\nPrinting the positions sampled, starting from the origin (0,0,0)\n";
    std::cout << "and using GAUSSIAN transition probability: "<<std::endl;
    hydrogen hydro(x,y,z,step100gaus, "100", "gaussian");

    std::vector<std::vector<double> > pos100(N), pos210(N);
    for(unsigned int i=0; i<N; ++i){
      hydro.metropolis();
      pos100[i]=hydro.get_position();
    }
    //used for tuning acceptance to 50%
    std::cout << "100 acceptance: " << hydro.acceptance() << std::endl;

    hydro.reset_metropolis(x,y,z,step210gaus, "210","gaussian");

    for(unsigned int i=0; i<N; ++i){
      hydro.metropolis();
      pos210[i]=hydro.get_position();
    }
    std::cout << "210 acceptance: " << hydro.acceptance() << std::endl;

    std::ofstream out100("100.dat");
    std::ofstream out210("210.dat");
    for(unsigned int i=0; i<N;++i){
        out100 << pos100[i][0] << " "
               << pos100[i][1] << " "
               << pos100[i][2] << "\n";
        out210 << pos210[i][0] << " "
               << pos210[i][1] << " "
               << pos210[i][2] << "\n";
    }
    out100.close();
    out210.close();
  }

  return 0;
}
