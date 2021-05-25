/*
Angela Riva
EXERCISE 07

07.1) Add the possibility of equilibration in MC simulation;
07.2) Include the computation of the radial distribution function;
*/

#include "montecarlo.h"

int main(int argc, char** argv) {
  std::string usage =
  "Usage: " + std::string(argv[0]) + " [option1] [option2] [option3] \n"+
  "This program performs a molecular dynamics simulation in NVT ensemble; \n"+
  "option1 (compulsory): \n" +
  "\t--equilibration      to start from 'config.0' configuration\n"+
  "\t--restart            to restart simulation from 'config.final'\n"+
  "\t                     named 'config.final'\n"+
  "option2 (compulsory): \n"+
  "Select the input file among the followings:"+
  "\t input.dat \n"+
  "\t input.solid \n"+
  "\t input.liquid \n"+
  "\t input.gas \n"+
  "option3 (NON compulsory):\n"+
  "\t--instant             to print istantaneous values\n";

  if(argc!=3 and argc!=4) std::cout << "Invalid option; \n \n" << usage << std::endl;

  if(argc==3) {
    if(std::string(argv[1])=="--restart") {
      MolDynMC MDMC(argv[2], "config.final");
      MDMC.Run();
    } else if(std::string(argv[1]) == "--equilibration"){
      MolDynMC MDMC(argv[2], "config.0");
      MDMC.Run();
    }
  } else if(argc==4 and std::string(argv[3])=="--instant") {
    if(std::string(argv[1])=="--restart") {
      MolDynMC MDMC(argv[2], "config.final", true);
      MDMC.Run();
    } else if(std::string(argv[1]) == "--equilibration"){
      MolDynMC MDMC(argv[2], "config.0", true);
      MDMC.Run();
    }
  } else std::cout << "Unrecognized option; \n\n " << usage << std::endl;


  return 0;
}
