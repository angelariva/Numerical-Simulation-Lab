#ifndef _TSP_
#define _TSP_

#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <algorithm>  // to use random_shuffle
#include <cmath>
#include <vector>
#include <utility>    // to use std::pair
#include <array>
#include "random.h"

using namespace std;
using city = std::pair<double, double>;   // a city is a pair of coordinates!
using path = std::vector<unsigned int>;   // the path of the salesmans

class SalesMan {
private:
	unsigned int Ncities;
	path Path;
	vector<city> Cities;
  Random* rnd;
  double pair_prob;
  double multi_prob;
  double inve_prob;
  double shift_prob;
public:

	// Constructors
	SalesMan();                        // In the main constructor:
	SalesMan(Random*,                  // POINTER to random object;
           string shape,             // shape of the world where the SalesMan travels;
           unsigned int ncities,     // total number of cities to be visited;
           double pair_prob=0.05,    // probability of pair permutation;
           double multi_prob=0.05,   // probability of permutation;
           double inve_prob=0.05,    // probability of inversion;
           double shift_prob=0.05);  // probability of shift.

	SalesMan(const SalesMan&);

	// Destructors
	~SalesMan();

	// Methods
	path GetPath();
	void PrintPath();
	void SetPath(path);

	vector<city> GetCities();          // In this method we get the coordinates of
                                     // the cities but in the order in which
                                     // the SalesMan visits them!
	void SetCities(vector<city>);

	bool Check();
  double AbsoluteLoss();             // Loss function defined as the sum of the
                                     // absolute distances between the cities
                                     // visited by the SalesMan.

	// Pbc
	int PbcPath(int);                  // Periodic Boundary Conditions for the
                                     // SalesMan's path.
	int PbcMutation(int);              // Will be used in applying the mutations.

  // Mutations
	void PairPermutation();
	void Shift();
	void MultiPermutation();
	void Inversion();

	// Simulated Annealing
	void SimulatedAnnealing(SalesMan&, 	// The starting SalesMan
													string, 					// the shape
													unsigned int, 		// # of temperatures to be considered
													unsigned int);		// # of steps
};

class Population {

private:
	unsigned int Npop;           // The # of individuals;
	unsigned int Ncities;        // The # of cities;
	vector<SalesMan> Pop;        // The population of SalesMen;
	vector<double> Losses;       // A vector containing the losses of the SalesMen;
  Random* rnd;                 // Pointer to random object;
	double cross_prob;           // Crossover probability.
public:

	// Constructors
                                     // The constructor needs:
	Population(Random*,                // Pointer to a random generator;
             string,                 // Shape of the world;
             unsigned int,           // # of salesmen;
             unsigned int,           // # of cities;
             double cross_prob=0.7); // crossover probability.
	// Destructors
	~Population();

	// Methods
	vector<SalesMan> GetPop();
	void SetPop(vector<path>, unsigned int);
	void PrintPop();
  void PrintLosses();

	SalesMan GetSalesMan(unsigned int);
	vector<double> GetLosses();
  double LossesAverage();

	// Crossover
	void OrderPop();
	unsigned int Selection();
	vector<path> Crossover();

	// Mutations
	void Mutations();

  // Evolution
  void Evolve(unsigned int, std::string shape); // Evolve the population.
};



#endif
