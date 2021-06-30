/*
Angela Riva

10.2)
  This program solves the Travelling SalesMan Problem using
  the mpi.h libraries.
*/

#include "TSP.h"

using namespace std;


int main (int argc, char *argv[]){

  std::string usage =
  "Usage: " + std::string(argv[0]) + " [Npop] [Ngen] [shape] \n" +
  "This program solves the Travelling SalesMan Problem using Genetic Algorithm \n"+
  "& using the mpi.h libraries for parallel computing. \n"+
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
    int Rank;
    int size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Status stat;

    double cross=0.7;
    Random rnd;
    rnd.SetRandom(Rank);
    srand(Rank + 1);  // we are using it in std::random_shuffle,
                      // we need it to be different for each run of the program!
    unsigned int ncit = 32;
    unsigned int nmigr = 10;
    unsigned int npop = atof(argv[1]);
    unsigned int ngen = atof(argv[2]);
    std::string shape = std::string(argv[3]);

  	ofstream cities("./results/cities_initial_"+shape+"_"+to_string(Rank)+".dat");

  	Population pop(&rnd, shape, npop, ncit, cross);
    vector<city> Cities(ncit);
    vector<double> xcoord(ncit);
    vector<double> ycoord(ncit);

    if(Rank==0) {
      for(unsigned int i=0; i<ncit; ++i) {
        xcoord[i] = pop.GetSalesMan(0).GetCities()[i].first;
        ycoord[i] = pop.GetSalesMan(0).GetCities()[i].second;
      }
    }

    MPI_Bcast(&xcoord.front(), ncit, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ycoord.front(), ncit, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for(unsigned int i=0; i<ncit; ++i){
			Cities[i].first = xcoord[i];
			Cities[i].second = ycoord[i];
    }

    pop.SetCities(Cities);

  	for(unsigned int i=0; i<ncit; ++i) {
  		cities << (pop.GetSalesMan(0).GetCities()[i]).first << "  "
             << (pop.GetSalesMan(0).GetCities()[i]).second << endl;
  	}
    cities << (pop.GetSalesMan(0).GetCities()[0]).first << "  "
           << (pop.GetSalesMan(0).GetCities()[0]).second << endl;

  	cities.close();

  	cities.open("./results/cities_final_"+shape+"_"+to_string(Rank)+".dat");

  	pop.OrderPop();

    if(Rank==0) {
      cout << "------------------------------" << endl;
      cout << ncit << " cities on a " << shape << ";" << endl;
      cout << "Population of " << npop << " salesmen;" << endl;
      cout << "Evolution of " << ngen << " generations. \n" << endl;
    }

    cout <<  "Rank_" << Rank << ") Best Initial loss value: "
         << pop.GetLosses()[npop - 1] << endl;

    pop.ParallelEvolve(Rank, ngen, nmigr, shape, &stat);

  	cout <<  "Rank_" << Rank << ") Best final loss value: "
         << pop.GetLosses()[npop-1] << endl;

  	for(unsigned int i=0; i<ncit; i++) {
  		cities << pop.GetSalesMan(npop-1).GetCities()[i].first << "  "
  			     << pop.GetSalesMan(npop-1).GetCities()[i].second << endl;
  	}

  	cities << pop.GetSalesMan(npop-1).GetCities()[0].first << "  "
  		     << pop.GetSalesMan(npop-1).GetCities()[0].second << endl;

    MPI_Finalize();
    cities.close();
  }
	return 0;
}
