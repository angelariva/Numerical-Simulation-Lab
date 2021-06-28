#include "TSP.h"

using namespace std;


int main (int argc, char *argv[]){

  std::string usage =
  "Usage: " + std::string(argv[0]) + " [Npop] [Ngen] [shape] \n" +
  "This program solves the Travelling SalesMan Problem using Genetic Algorithm. \n"+
  "The (mandatory) options to be inserted from command line are: \n"+
  "\t [Npop]:  # of individuals in the first generation; \n"+
  "\t [Ngen]:  # of generations; \n"+
  "\t [shape]: shape of the world. The only two possibilities are: \n"+
  "\t\t  - circumference; \n "+
  "\t\t  - square.";

  if(argc!=4) {
    std::cout << "Invalid Option." << std::endl;
    std::cout << usage << std::endl;
  } else {
    double cross=0.7;
    Random rnd;
    rnd.SetRandom("Primes", "seed.in");
    unsigned int ncit = 32;
    unsigned int npop = atof(argv[1]);
    unsigned int ngen = atof(argv[2]);
    std::string shape = std::string(argv[3]);

  	ofstream cities("./results/cities_"+shape+".dat");
  //	ofstream losses("./results/cost_"+shape+".dat");
  //	ofstream losses_ave("./results/cost_ave_"+shape+".dat");

  	Population pop(&rnd, shape, npop, ncit, cross);
  	vector<path> appop(2);
  	vector<path> popul(npop);

  	for(unsigned int i=0; i<ncit; ++i) {
  		cities << (pop.GetSalesMan(0).GetCities()[i]).first << "  "
             << (pop.GetSalesMan(0).GetCities()[i]).second << endl;
  	}
  	cities.close();

  	cities.open("./results/cities_final_"+shape+".dat");

  	pop.OrderPop();

    cout << "------------------------------" << endl;
    cout << ncit << " cities on a " << shape << ";" << endl;
    cout << "Population of " << npop << " salesmen;" << endl;
    cout << "Evolution of " << ngen << " generations. \n" << endl;
  	cout << "Best Initial loss value: " << pop.GetLosses()[npop - 1] << endl;

    pop.Evolve(ngen, shape);

  	cout << "Best final loss value: " << pop.GetLosses()[npop-1] << endl;

  	for(unsigned int i=0; i<ncit; i++) {
  		cities << pop.GetSalesMan(npop-1).GetCities()[i].first << "  "
  			     << pop.GetSalesMan(npop-1).GetCities()[i].second << endl;
  	}

  	cities << pop.GetSalesMan(npop-1).GetCities()[0].first << "  "
  		     << pop.GetSalesMan(npop-1).GetCities()[0].second << endl;

  	cities.close();
//  	losses.close();
//  	losses_ave.close();
  }
	return 0;
}
