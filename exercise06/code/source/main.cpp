/*
Angela Riva
EXERCISE 06
this program can be executed with different options.
It simulates the behavior of an Ising model. 
*/

#include "ising.h"

int main(int argc, char** argv) {
    std::string usage =
    "Usage: " + std::string(argv[0]) + " [options] \n"+
    "Where options (compulsory) should be: \n" +
    "\t--first             to start from random configuration\n"+
    "\t--restart           to restart simulation from configuration file\n"+
    "\t                    named 'config.final'\n"+
    "\t--equilibration     to start from random config and print \n"+
    "\t                    the instantaneous values of the equilibration";

    if(argc!=2){
        std::cout << "Program for simulating 1D Ising Model\n\n" <<usage;
        return 0;
    }

    if(std::string(argv[1])=="--restart") {
        Ising1D ising("config.final");
        ising.Run();
    }
    else if(std::string(argv[1])=="--equilibration") {
        Ising1D ising;
        ising.Run(true);
    }
    else if(std::string(argv[1])=="--first") {
        Ising1D ising;
        ising.Run();
    }
    else
        std::cout << "\nUnrecognized option '" << argv[1] << "'\n" <<
            std::endl <<usage;
    return 0;
}
