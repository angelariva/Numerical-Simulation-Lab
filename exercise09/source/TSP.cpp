#include "TSP.h"


///////////////////////////////////////
//       SalesMan                    //
///////////////////////////////////////

SalesMan::SalesMan() { }

SalesMan::SalesMan(Random* r,
                   string shape,
                   unsigned int ncities,
                   double pair_prob,
                   double multi_prob,
                   double inve_prob,
                   double shift_prob) : Cities(ncities),
                                        pair_prob(pair_prob),
                                        multi_prob(multi_prob),
                                        inve_prob(inve_prob),
                                        shift_prob(shift_prob) {
  rnd = r;
	Ncities = ncities;
	for(unsigned int i=1; i<=Ncities; ++i) {
		Path.push_back(i);
	}

	random_shuffle(Path.begin() + 1, Path.end());

 	if(shape=="circumference"){
		double theta;
		for(unsigned int i=0; i<Ncities; ++i){
			theta = rnd->Rannyu(0, 2*M_PI);
			Cities[i].first = cos(theta);
			Cities[i].second = sin(theta);
		}
	} else if(shape=="square") {
		for(unsigned int i=0; i<Ncities; ++i){
			Cities[i].first = 2 * rnd->Rannyu() - 1.;
			Cities[i].second = 2 * rnd->Rannyu() - 1.;
		}
	}

}

SalesMan::~SalesMan() { }

path SalesMan::GetPath() {
	return Path;
}

void SalesMan::PrintPath() {
	cout << "[";
	for(unsigned int i=0; i<Ncities; ++i) {
		cout << Path[i];
		if(i < Ncities - 1) cout << ",";
	}
	cout << "]";
}

void SalesMan::SetPath(path p) {
	Path = p;
}

vector<city> SalesMan::GetCities() {
	vector<city> cities(Ncities);
	for(unsigned int i=0; i<Ncities; ++i){
		cities[i] = Cities[Path[i]-1];
	}
	return cities;
}

void SalesMan::SetCities(vector<city> cities) {
	Cities = cities;
}

bool SalesMan::Check() {

	bool check=true;
	int count;
	if(Path[0] != 1) check = false;
	else{
		for(unsigned int i=1; i<=Ncities; ++i){
			count = 0;
			for(unsigned int j=0; j<Ncities; j++) if(Path[j] == i) count++;
			if(count > 1) { check = false; break;}
		}
	}
	return check;
}

double SalesMan::AbsoluteLoss() {

	double loss = 0;
	for(unsigned int i=0; i<Ncities; ++i){
		loss += sqrt(pow(Cities[Path[i]-1].first - Cities[Path[PbcPath(i+1)]-1].first, 2)
			 + pow(Cities[Path[i]-1].second - Cities[Path[PbcPath(i+1)]-1].second, 2));
	}

	return loss;

}

int SalesMan::PbcPath(int i) {
	return i % Ncities;
}

int SalesMan::PbcMutation(int i) {
	if(abs(i)<Ncities) return i;
	else			 return i % Ncities + 1;
}

void SalesMan::PairPermutation() {
	unsigned int i1, i2;
	if(rnd->Rannyu() < pair_prob) {
		i1 = rnd->RandInt(1, Ncities - 1);
		do{
			i2 = rnd->RandInt(1, Ncities - 1);
		} while(i2 == i1);
		swap(Path[i1], Path[i2]);
	}
}

void SalesMan::Shift() {
	unsigned int i1, m, n;
	if(rnd->Rannyu() < shift_prob){
		i1 = rnd->RandInt(1, Ncities - 1);
		m = rnd->RandInt(1, Ncities - 1);

		// n goes from 1 to (m-1). At n = m
		// we go back to the initial SalesMan
		n = rnd->RandInt(1, m);
		for(unsigned int i=0; i<n; ++i){
			// This for represent a single shift
			for(unsigned int j=0; j<m-1; j++){
				swap(Path[PbcMutation(i1 + j)], Path[PbcMutation(i1 + m - 1)]);
			}
		}


	}
}

void SalesMan::MultiPermutation() {
	unsigned int m, i1, i2;
	if(rnd->Rannyu() < multi_prob) {
		m = rnd->RandInt(1, (Ncities + 1) / 2 - 1); // 1 <= m < Ncities/2
		i1 = rnd->RandInt(1, Ncities - 1);
		i2 = PbcMutation(rnd->RandInt(i1 + m, Ncities + i1 - m - 1));
		for(unsigned int i=0; i<m; ++i) {
			swap(Path[PbcMutation(i1+i)], Path[PbcMutation(i2+i)]);
		}
	}
}

void SalesMan::Inversion() {
	int i1, m;
	if(rnd->Rannyu() < inve_prob) {
		i1 = rnd->RandInt(1, Ncities - 1);
		m = rnd->RandInt(1, Ncities);
		for(unsigned int i=0; i<abs((int)m/2); ++i) {
			swap(Path[PbcMutation(i1 + i)], Path[PbcMutation(i1 + (m - 1) - i)]);
		}
	}
}

///////////////////////////////////////
//       Population                  //
///////////////////////////////////////

Population::Population(Random* r,
                       string shape,
                       unsigned int npop,
                       unsigned int ncities,
                       double cross) : Losses(npop, 0) {
  rnd  = r;
  cross_prob = cross;
	Npop = npop;
	Ncities = ncities;

	for(unsigned int i=0; i<Npop; ++i) {
		Pop.push_back(SalesMan(r, shape, Ncities));
	}

	vector<city> cities(Ncities);

  if(shape=="circumference") {
	double theta;
	for(unsigned int i=0; i<Ncities; ++i) {
			theta = rnd->Rannyu(0, 2*M_PI);
			cities[i].first = cos(theta);
			cities[i].second = sin(theta);
		}
		for(unsigned int i=0; i<Npop; ++i) {
			Pop[i].SetCities(cities);
		}
	} else if(shape=="square") {
		for(unsigned int i=0; i<Ncities; ++i) {
			cities[i].first = 2 * rnd->Rannyu() - 1.;
			cities[i].second = 2 * rnd->Rannyu() - 1.;
		}
		for(unsigned int i=0; i<Npop; ++i) {
			Pop[i].SetCities(cities);
		}
	}
}

Population::~Population() { }

SalesMan Population::GetSalesMan(unsigned int i) {
	return Pop[i];
}

vector<SalesMan> Population::GetPop() {
	return Pop;
}

void Population::PrintPop(){
	cout << "------------------------------------------------" <<endl;
	for(unsigned int i=0; i<Npop; ++i) {
		Pop[i].PrintPath();
		cout << "\n AbsoluteLoss = " << Pop[i].AbsoluteLoss() << endl;
	}
	cout << "------------------------------------------------" <<endl;
}

void Population::SetPop(vector<path> pop, unsigned int n){
	for(unsigned int i=0; i<n; ++i) {
		Pop[i].SetPath(pop[i]);
	}
}

vector<double> Population::GetLosses(){
	return Losses;
}

double Population::LossesAverage() {
  double ave = 0.;
	for(unsigned int k = Npop-1; k>Npop/2; k--) {
		ave += Losses[k];
	}
	ave /= (Npop/2 - 1);
	return ave;
}

struct Mileage {
	SalesMan Chromo;
	double Miles;
};

void Population::OrderPop(){

	vector<Mileage> mileage(Npop);
	for(unsigned int i=0; i<Npop; ++i){
		mileage[i].Chromo = Pop[i];
		mileage[i].Miles = Pop[i].AbsoluteLoss();
	}

	// We order the population on a fitness basis:
	// the higher AbsoluteLoss, the lower the rank
	sort(mileage.begin(), mileage.end(),
		 [](Mileage const &a, Mileage const &b) { return a.Miles > b.Miles; });

	// We write the population (and the respective mileage) in a sorted fashion
	 for(unsigned int i=0; i<Npop; ++i) {
        Pop[i] = mileage[i].Chromo;
        Losses[i] = mileage[i].Miles;
    }
}

void Population::PrintLosses() {
  for(auto& el : Losses) cout << "Losses: " << el << endl;
  cout << "----------------------------------------- " << endl;
}

unsigned int Population::Selection() {
	double p = 1./5;
	double r = rnd->Rannyu();
	return int(Npop * pow(r, p));
}

vector<path> Population::Crossover() {
	vector<path> offspring(2);
	path offspring1, offspring2;
	unsigned int i1, i2, cut;
  unsigned int start1=1, start2 = 1;
  unsigned int elem;

  i1 = Selection();
  do {
    i2 = Selection();
  } while(i2 == i1);

  if(rnd->Rannyu() < cross_prob) {
    cut = rnd->RandInt(1, Ncities - 1);
    // At this position we cut the SalesMans

    offspring1 = Pop[i1].GetPath();
    offspring2 = Pop[i2].GetPath();

    for(unsigned int i=cut; i<Ncities; ++i) {
      // We use this cycle to write offspring1 properly
      for(unsigned int j=start1; j<Ncities; ++j) {
        elem = Pop[i2].GetPath()[j];
        if(find(offspring1.begin(), offspring1.begin() + cut, elem) == offspring1.begin() + cut) {
          // if offspring2[j] is not in offspring1 from 0 to cut-1
          offspring1[i] = elem;
          start1++;
          break;
        } else { start1++;}
      }
      // We use this cycle to write offspring2 properly
      for(unsigned int j=start2; j<Ncities; ++j){
        elem = Pop[i1].GetPath()[j];
        if(find(offspring2.begin(), offspring2.begin() + cut, elem) == offspring2.begin() + cut) {
          // if offspring1[j] is not in offspring2 from 0 to cut-1
          offspring2[i] = elem;
          start2++;
          break;
        } else {start2++;}
      }
    }
    offspring[0] = offspring1;
    offspring[1] = offspring2;
  } else {
    // If there is no crossover we just copy
    // the selected SalesMans
    offspring[0] = Pop[i1].GetPath();
    offspring[1] = Pop[i2].GetPath();
  }
  return offspring;
}

void Population::Evolve(unsigned int gen, std::string shape) {
  ofstream loss("./results/loss_"+shape+".dat");
  ofstream lossave("./results/loss_ave_"+shape+".dat");
  vector<path> offsprings(2);
  vector<path> pop(Npop);

  for(unsigned int j=0; j<gen; ++j) {
		for(unsigned int i=0; i<Npop/2; ++i) {
        offsprings = (*this).Crossover();
		    pop[2*i] = offsprings[0];
		    pop[2*i + 1] = offsprings[1];
		}

		(*this).SetPop(pop, Npop);
		(*this).Mutations();
		(*this).OrderPop();
//    (*this).PrintLosses();
//    (*this).PrintPop();
    loss << j+1 << " " << GetLosses()[Npop - 1] << endl;
		lossave << j+1 << " " << LossesAverage() << endl;

	}
  loss.close();
  lossave.close();
}

void Population::Mutations() {
	for(unsigned int i=0; i<Npop; ++i) {
		Pop[i].PairPermutation();
		Pop[i].Shift();
		Pop[i].MultiPermutation();
		Pop[i].Inversion();
	}
}
