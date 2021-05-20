/*
Angela Riva
EXERCISE 04
In this program we exploit the class MolecularDynamics.
Depending on the input option, at the moment of the execution of the program,
a different simulation is run:
- the equilibration simulation starting from the zero configuration;
- the actual simulation starting from old configurations.
This is done thanks to the two diffenrent constructors implemented in the class.
*/
#include "MolDyn_NVE.h"

using namespace std;

void usage() {
  cout << endl;
  cout << "This program performs a molecular dynamics simulation in NVE ensemble;" << endl;
  cout << "Usage : " << "./main [inputfile] [option2] " << endl;
  cout << "option1: " << endl;
  cout << "   --equilibration     " << ": starts from the configuration in config.0"  << endl;
  cout << "   --restart           " << ": starts from configuration in old.0 " << endl;
  cout << "                               and config.final" << endl << endl;
  cout << "option2: select the input file among the followings: " << endl;
  cout << "   input.dat"  << endl;
  cout << "   input.solid"  << endl;
  cout << "   input.liquid"  << endl;
  cout << "   input.gas"  << endl;
}

int main (int argc, char *argv[]) {
  // argc: # of strings appearing on command line
  // argv: array of strings, argv[0] = ./exec_name

  if(argc!=3) usage();

  if(std::string(argv[1])=="--restart"){
    MolecularDynamics molDin(argv[2], "config.final", "old.0");
    molDin.RunSimulation();
  } else if(std::string(argv[1])=="--equilibration"){
    MolecularDynamics molDin(argv[2], "config.0");
    molDin.RunSimulation();
  } else usage();

  return 0;
}
