/*
Angela Riva

10.1)
  This program solves the Travelling SalesMan Problem using
  Simulated Annealing.
*/

#include "TSP.h"

using namespace std;


int main (int argc, char *argv[]) {

  std::string usage =
  "Usage: " + std::string(argv[0]) + " [Ntemp] [Nstep] [shape] \n" +
  "This program solves the Travelling SalesMan Problem using Genetic Algorithm. \n"+
  "The (mandatory) options to be inserted from command line are: \n"+
  "\t [Ntemp]:  # of temperatures to be considered; \n"+
  "\t [Nstep]:  # of steps to be made for each temperature; \n"+
  "\t [shape]: shape of the world. The only two possibilities are: \n"+
  "\t\t  - circumference; \n "+
  "\t\t  - square.";

  if(argc!=4) {
    std::cout << "Invalid Option." << std::endl;
    std::cout << usage << std::endl;
  } else {
    Random rnd;
    rnd.SetRandom("Primes", "seed.in");
    unsigned int ngen = 200;
    unsigned int npop = 1000;
    unsigned int ncit = 32;
    unsigned int ntemp = atof(argv[1]);
    unsigned int nstep = atof(argv[2]);
    std::string shape = std::string(argv[3]);

  	ofstream cities("./results/cities_"+shape+".dat");
    ofstream res("./results/results_"+shape+".txt");

    Population pop(&rnd, shape, npop, ncit, 0.7);
    vector<city> coord = (pop.GetSalesMan(0)).GetCities();
    pop.OrderPop();
    cout << "------------------------------" << endl;
    cout << "Genetic Algortithm" << endl;
    cout << ncit << " cities on a " << shape << ";" << endl;
    cout << "Population of " << npop << " salesmen;" << endl;
    cout << "Evolution of " << ngen << " generations. \n" << endl;
    cout << "Best Initial loss value: " << pop.GetLosses()[npop - 1] << endl;

    res << "------------------------------" << endl;
    res << "Genetic Algortithm" << endl;
    res << ncit << " cities on a " << shape << ";" << endl;
    res << "Population of " << npop << " salesmen;" << endl;
    res << "Evolution of " << ngen << " generations. \n" << endl;
    res << "Best Initial loss value: " << pop.GetLosses()[npop - 1] << endl;

    pop.Evolve(ngen, shape);
    cout << "Best final loss value: " << pop.GetLosses()[npop-1] << endl;
    res << "Best final loss value: " << pop.GetLosses()[npop-1] << endl;

    SalesMan Initial(&rnd, shape, ncit);
    Initial.SetCities(coord);
    SalesMan Final = Initial;
    for(unsigned int i=0; i<ncit; ++i) {
      cities << Initial.GetCities()[i].first << "  "
             << Initial.GetCities()[i].second << endl;
  	 }

    cities.close();
  	cities.open("./results/cities_final_"+shape+".dat");

    cout << "------------------------------" << endl;
    cout << "Simulated Annealing" << endl;
    cout << ncit << " cities on a " << shape << ";" << endl;
    cout << ntemp << " temperatures to be considered;" << endl;
    cout << nstep << " steps to be made for each temperature. \n"<< endl;
  	cout << "Initial loss value: " << Initial.AbsoluteLoss() << endl;

    res << "------------------------------" << endl;
    res << "Simulated Annealing" << endl;
    res << ncit << " cities on a " << shape << ";" << endl;
    res << ntemp << " temperatures to be considered;" << endl;
    res << nstep << " steps to be made for each temperature. \n"<< endl;
  	res << "Initial loss value: " << Initial.AbsoluteLoss() << endl;

    Final.SimulatedAnnealing(Initial, shape, ntemp, nstep);

    cout << "Final loss value: " << Initial.AbsoluteLoss() << endl;
    res << "Final loss value: " << Initial.AbsoluteLoss() << endl;

    for(unsigned int i=0; i<ncit; ++i) {
		    cities << Initial.GetCities()[i].first << "  "
               << Initial.GetCities()[i].second << endl;
	  }

    cities << Initial.GetCities()[0].first << "  "
           << Initial.GetCities()[0].second << endl;
  	cities.close();

  }
	return 0;
}
