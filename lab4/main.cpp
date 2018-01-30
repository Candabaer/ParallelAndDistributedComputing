#define _SCL_SECURE_NO_WARNINGS

// g++ -c main.cpp -pthread -std=c++11 -fopenmp -O3
// g++ main.o -pthread -std=c++11 -fopenmp -O3
// ./a.out

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <chrono>
#include <iterator>
#include <set>
#include <omp.h>
#include <thread>
#include <limits.h>
#include <list>

using namespace std;

int NT = 4;

bool abortEvol;

void timer(int min) {
	std::this_thread::sleep_for(std::chrono::minutes(min));
	abortEvol = true;
}

string dts(double d) {
	std::ostringstream strs;
	strs << d;
	string s = strs.str();
	return s;
}

//Node structure
class City {
public:
	City() {}
	City(int ID, string name, string plz, double longi, double lat) {
		this->ID = ID;
		this->name = name;
		this->plz = plz;
		this->longi = longi;
		this->lat = lat;
		this->visited = false;
	}

	City(City* a) {
		this->ID = a->ID;
		this->name = a->name;
		this->plz = a->plz;
		this->longi = a->longi;
		this->lat = a->lat;
		this->visited = a->visited;
	}

	~City() {}

	int distanceToCity(City* city) {
		double RRR = 6378.388;
		double q1 = cos(this->longi - city->longi);
		double q2 = cos(this->lat - city->lat);
		double q3 = cos(this->lat + city->lat);
		return (int)(RRR * acos(0.5*((1.0 + q1)*q2 - (1.0 - q1)*q3)) + 1.0);
	}

	string getName() {
		return this->name;
	}
	int getId() {
		return this->ID;
	}

	void setID(int id) {
		this->ID = id;
	}

	double getLat() {
		return this->lat;
	}

	double getLongi() {
		return this->longi;
	}

	void printCity() {
		cout << "ID: " << this->ID << ", name: " << this->name << ", plz: " << this->plz << ", longi: " << this->longi << ", lat: " << this->lat << endl;
	}
private:
	int ID;
	string name;
	string plz;
	double longi;
	double lat;
	bool visited;
};

//Holds every existing Node
class Map {
public:
	vector<City*> nodes;
	int size;
	int** distMatrix;

	Map(string fname) {
		ifstream is(fname);
		if (!is) {
			cerr << "Error, could not open file: " << fname << endl;
		}
		char c[256];
		is.getline(c, 256);
		while (is.good()) {
			std::stringstream ss;
			int id;
			std::string name;
			std::string plz;
			double lonitude;
			double latidude;
			ss << c;
			ss >> id >> plz >> lonitude >> latidude >> name;
			City* newNode = new City(id, name, plz, lonitude, latidude);
			this->nodes.push_back(newNode);
			is.getline(c, 256);
		}
		this->size = nodes.size();
		this->calcDistMatrix();
		cout << "Added and Calculated: " << this->size << " Citys." << endl;	
		// this->printCitys();
	}

	~Map(){}
	
	void printCitys() {
		for (int z = 0; z < this->size; z++) {
			this->nodes[z]->printCity();
		}
	}

private: 
	void calcDistMatrix() {
		int** distanceMatrix = new int*[this->size];
		for (int z = 0; z < this->size; z++) {
			distanceMatrix[z] = new int[this->size];
		}
#pragma omp parallel for
		for (int i = 0; i < this->size; i++) {
			// distanceMatrix[i] = new int[this->size];
			for (int ii = 0; ii <this->size ; ii++) {
				if (i == ii) {
					distanceMatrix[i][ii] = INT_MAX;
					continue;
				}
				double RRR = 6378.388;
				double q1 = cos(this->nodes[i]->getLongi() - this->nodes[ii]->getLongi());
				double q2 = cos(this->nodes[i]->getLat() - this->nodes[ii]->getLat());
				double q3 = cos(this->nodes[i]->getLat() + this->nodes[ii]->getLat());
				distanceMatrix[i][ii] = (int)(RRR * acos(0.5*((1.0 + q1)*q2 - (1.0 - q1)*q3)) + 1.0);
			}
		}
		this->distMatrix = distanceMatrix;
	}
};

//Holds all Nodes in a given Order as a Tour.
class Tour {
private:
	int* cities;
	int size;
	int distance;
	double fitness;
	bool changed;

	void calcFitness() {
		this->fitness = 10000000000*(1 / (double)this->distance);
	}
public:
	Tour() {
		this->size = 0;
		this->cities = NULL;
		this->fitness = 0;
		this->distance = 0;
		this->changed = true;
	}

	Tour(int size) {
		this->size = size;
		this->cities = new int[this->size];
		this->fitness = 0;
		this->distance = 0;
		this->changed = true;
	}

	Tour(Tour* a) {
		this->changed = a->changed;
		this->size = a->size;
		this->distance = a->distance;
		this->fitness = a->fitness;
		this->cities = new int[this->size];
		for (int z = 0; z < this->size;z++) {
			this->cities[z] = int(a->cities[z]);
		}
	}

	~Tour() {
		 delete[] cities;
	}

	double getFitness() {
		return this->fitness;
	}

	double CalcFitnessAndDist(Map* map) {
		int dist = 0.;
		#pragma omp parallel for reduction (+:dist)
		for (int z = 0; z < this->size; z++) {
			int iFrom = this->cities[z];
			int iTo = z + 1;
			if (iTo < this->size) {
				iTo = this->cities[iTo];
			}
			else {
				iTo = this->cities[0];
			}
			//if (iFrom == iTo) {
			//	cout << "aaargh" << endl;
			//}
			dist = dist+map->distMatrix[iFrom][iTo];
		}
		this->distance = dist;
		this->calcFitness();
		return this->fitness;
	}

	int getDistance() {
		return this->distance;
	}

	bool containsCity(int id) {
		for (int z = 0; z < this->size; z++) {
			if (id == this->cities[z]) {
				return true;
			}
		}
		return false;
	}

	//Each Id in City* represents its place in the map->nodes[id]. Citys only exist in the Map!.
	void generateTour(Map* map) {
		for (int z = 0; z < this->size; z++) {
			this->cities[z] = int(z); // new City(map->nodes[z]);
		}
		std::random_shuffle(&this->cities[0],&this->cities[this->size]);
		
	}

	int* getCities() {
		return this->cities;
	}

	int getSize() {
		return this->size;
	}

	int getCity(int pos) {
		return this->cities[pos];
	}

	//searches for City with ID in range [r1, r2], if r2 is not specified, this->size is taken.
	int getCityByID(int ID, int r1 = 0, int r2 = 0) {
		if (r2 == 0) {
			r2 = this->size - 1;
		}
		for(int z=r1; z <= r2; z++)
			if (this->cities[z] == ID) {
				return this->cities[z];
			}
		std::cout << "No City found with ID: " << ID << ", in Range [" << r1 << "," << r2 << "]" << endl;
		return NULL;
	}

	void setCity(int pos, int ID) {
		this->cities[pos] = ID;
	}

	string getTourString(Map* map) {
		string str = "";
		for (int z = 0; z < this->size; z++) {
			str += std::to_string(map->nodes[this->cities[z]]->getId()) + " ";
		}
		return str;
	}

	//Returns the Index of the Citys place in MAP, not the actual City ID.
	void print() {
		std::cout << "Fitness: " << this->fitness << std::endl;
		for (int z = 0; z < this->size; z++) {
			cout << this->cities[z] << ", ";
		}
		cout << endl;
	}
	void hasError() {
		for (int z = 0; z < this->size; z++) {
			int id = this->cities[z];
			for (int s = 0; s < this->size; s++) {
				if (z == s) {
					continue;
				}
				if (id == this->cities[s]) {
					cout << "Error Found! " << z << " equals: "<< s <<endl;
				}
			}
		}
	}
	bool getChanged() {
		return this->changed;
	}

	void setChanged(bool to) {
		this->changed = to;
	}
};

class Population {
private:
	Tour** tours;
	Map* map;
	int size;
	int tourSize;
	double avgFitness;
	
public:
	// creates uninitialized, empty population
	Population(int size, int tourSize, Map* map) {
		this->size = size;
		this->tourSize = tourSize;
		this->tours = new Tour*[size];
		this->avgFitness = 0;
		this->map = map;
	}

	//Creates rand() initialised population with Tours out of Tour
	Population(Map* map, int size, int tourSize) {
		this->avgFitness = 0;
		this->size = size;
		this->tourSize = tourSize;
		this->tours = new Tour*[size];
		this->map = map;
		for (int z = 0; z < size; z++) {
				Tour* a = new Tour(tourSize);
				a->generateTour(map);
				this->tours[z] = a;
		}
		
	}

	~Population() {
		for (int z = 0; z < this->size; z++) {
			delete this->tours[z];
		}
		delete[] tours;
	}

	void saveTour(int index, Tour* tour) {
		this->tours[index] = tour;
	}

	Tour* getFittest() {		
		Tour* fittest = this->tours[0];
		for (int z = 1; z < this->size; z++) {
			if (fittest->getFitness() < this->tours[z]->getFitness()) {
				#pragma omp critical
				{
					fittest = this->tours[z];
				}
			}
		}
		return fittest;
	}

	void printAll() {
		for (int z = 0; z < this->size; z++) {
			cout << "ID: " << z << ", ";
			this->tours[z]->print();
		}
	}

	int getSize() {
		return this->size;
	}

	Tour* getTour(int index) {
		return this->tours[index];
	}

	double calcAverageFitness() {
		double fit = 0.;
		//#pragma omp parallel for reduction (+:fit)
		for (int z = 0; z < this->size; z++) {
			Tour* t = this->tours[z];
			if (t->getChanged() == true) {
				fit = fit + t->CalcFitnessAndDist(this->map);
			}
			else {
				fit = fit + t->getFitness();
			}			
		}
		this->avgFitness = fit / (double)this->size;
		return this->avgFitness;
	}

	void checkForErrorTours() {
		bool a = false;
		for (int z = 0; z < this->size; z++) {
			this->tours[z]->hasError();
		}
	}

	void setNull() {
	for (int z = 0; z < this->size; z++) {
		this->tours[z] = NULL;
	}
	this->map = NULL;
}

	double getAvgFitness() {
	return this->avgFitness;
}
};

class Statistics {
private:
	ofstream logFile;
public:
	Statistics() : logFile("data_output.csv", ofstream::trunc) {

	}
	void logData(int gen, double fitness, Tour* tour) {
		this->logFile << gen << ", " << fitness <<", , "<<gen << ", " << tour->getFitness() << endl;
	}

	void logTopCandidate(int gen, Tour* tour, Map* map) {
		this->logFile << endl << gen << ", " << tour->getFitness() << ", " << tour->getDistance() << ", " << tour->getTourString(map) << endl;
		ofstream tmp("candidateFromTo");
		for (int z = 0; z < tour->getSize(); z++) {
			if (z == tour->getSize()) {
				tmp << tour->getCity(z) << " " << tour->getCity(0) << endl;
			}
			else {
				tmp << tour->getCity(z) << " " << tour->getCity(z+1) << endl;
			}
		}

	}

	void logBasics(int popSize, int gens, double mut, double cross, double die) {
		this->logFile << "PopSize: " << popSize << ",MutationRate: " <<dts(mut) << ",Crossover: " << dts(cross) << ",dying: " << dts(die) <<",Gens: " << gens << endl;
	}

	void logTime(int seks) {
		this->logFile << "Runtime in Minutes.: " << seks << endl;
	}
};

class genAlgo {
private:
	double mutationRate = 0.05;
	double crossOverRate = 0.10;
	double dying = 0.05;
	int tournamentSize = 2;
	bool elitism = false;
	Population* pop;
	//Map* map = new Map("test_input.txt");
	Map* map = new Map("zips.txt");
	int popSize = 200;
	int gens = 0;
	Statistics* stats;
	std::pair<int, Tour*> bestTour;

	void mutate(Population* newPop, int toursize) {
		#pragma omp parallel for
		for (int i = 0; i < this->popSize; i++) {
			double rng = ((double)rand() / (RAND_MAX));
			if (rng < this->mutationRate) {				
				int pos1 = rand() % toursize;
				int pos2 = rand() % toursize;
				while (pos1 == pos2) {
					pos2 = rand() % toursize;
				}				
				int* cities = newPop->getTour(i)->getCities();
				#pragma omp critical
				{
				int val1 = cities[pos1];
				int val2 = cities[pos2];
				cities[pos1] = val2;
				cities[pos2] = val1;
				}
				newPop->getTour(i)->setChanged(true);
			}
		}
		//tour->print();
		//int pos1 = rand() % toursize;
		//int pos2 = rand() % toursize;	
		//if (pos1 > pos2) {
		//	int swt = pos2;
		//	pos2 = pos1;
		//	pos1 = swt;
		//}
		//int span = pos2 - pos1 + 1;
		//int* cities = tour->getCities();
		//int* arr = new int[span];
		//copy(&cities[pos1], &cities[pos2 + 1], &arr[0]);
		//std::reverse_copy(&arr[0], &arr[span], &cities[pos1]);
		//delete arr;
	}

	void selectCrossOver(Population* newPop, int toursize) {
		int amount = this->popSize * this->crossOverRate;
		for (int z = 0; z < amount; z++) {
			int a = rand() % this->popSize;
			int b;
			do {
				b = rand() % this->popSize;
			} while (a == b);
			crossOver(newPop, a, b, toursize);
		}
	}

	void crossOver(Population* newPop, int idA, int idB, int toursize) {
		Tour* a = newPop->getTour(idA);
		Tour* b = newPop->getTour(idB);
		//cout << "A: " << endl;
		//a->print();
		//cout << "B: " << endl;
		//b->print();
		Tour* childAMidB = new Tour(a);
		Tour* childBMidA = new Tour(b);
		int pos1 = rand() % toursize;
		int pos2;
		do {
			pos2 = rand() % toursize;
		} while (pos2 == pos1);
		int* citiesA = a->getCities();
		int* citiesB = b->getCities();
		int* citiesCAMidB = childAMidB->getCities();
		int* citiesCBMidA = childBMidA->getCities();
		if (pos1 > pos2) {
			int swt = pos2;
			pos2 = pos1;
			pos1 = swt;
		}		
		// find Mapping
		vector<int> vecA;
		vector<int> vecB;
		for (int z = pos1; z <= pos2; z++) {
			vecA.push_back(citiesA[z]);
			vecB.push_back(citiesB[z]);
			citiesCAMidB[z] = int(citiesB[z]);
			citiesCBMidA[z] = int(citiesA[z]);
		}
		for (int z = 0; z < vecA.size(); z++) {
			std::vector<int>::iterator it = find(vecB.begin(), vecB.end(), vecA[z]);
			if (it != vecB.end()) {	//element gefunden!
				*it = vecB[z];
				vecA.erase(vecA.begin() + z);
				vecB.erase(vecB.begin() + z);
				z--;
			}
		}

		//cout << "Child before Apllying Map: " << endl;
		//child->print();
		/*cout << "Crossover in: " << pos1 << " - " << pos2 << ". Mapping: " << endl;
		for (int z = 0; z < vecA.size(); z++) {
			cout << vecA[z] << " -> " << vecB[z] << endl;
		}*/
		// #pragma omp parallel for
		for (int z = 0; z < pos1; z++) {
			int cIDAMidB = citiesA[z];
			int cIDBMidA = citiesB[z];
			for (int m = 0; m < vecA.size(); m++) {
				if (cIDAMidB == vecA[m]) {
						citiesCAMidB[z] = int(b->getCityByID(vecB[m],pos1,pos2));
				   //break;
				}
				else if (cIDAMidB == vecB[m]) {
					citiesCAMidB[z] = int(a->getCityByID(vecA[m],pos1,pos2));
				   //break;
				}
				if (cIDBMidA == vecA[m]) {
					citiesCBMidA[z] = int(b->getCityByID(vecB[m], pos1, pos2));
				}
				else if (cIDBMidA == vecB[m]) {
					citiesCBMidA[z] = int(a->getCityByID(vecA[m], pos1, pos2));
				}
			}
		}
		// #pragma omp parallel for
		for (int z = pos2+1; z < toursize; z++) {
			int cIDAMidB = citiesA[z];
			int cIDBMidA = citiesB[z];
			for (int m = 0; m < vecA.size(); m++) {
				if (cIDAMidB == vecA[m]) {
					citiesCAMidB[z] = int(b->getCityByID(vecB[m], pos1, pos2));
					// break;
				}
				else if (cIDAMidB == vecB[m]) {
					citiesCAMidB[z] = int(a->getCityByID(vecA[m], pos1, pos2));
					// break;
				}
				if (cIDBMidA == vecA[m]) {
					citiesCBMidA[z] = int(b->getCityByID(vecB[m], pos1, pos2));
				}
				else if (cIDBMidA == vecB[m]) {
					citiesCBMidA[z] = int(a->getCityByID(vecA[m], pos1, pos2));
				}
			}
		}
		// cout << "After mapping: AMidB " << endl;
		// childAMidB->print();
		//cout << "After mapping: BMidA " << endl;
		// childBMidA->print();
		delete a;
		delete b;
		newPop->saveTour(idA, childAMidB);
		newPop->saveTour(idB, childBMidA);
		//return childAMidB;
	}
	
	void select(Population* newPop, int toursize) {
		#pragma omp parallel for
		for (int i = 0; i < this->popSize; i++) {
			double rngJesus = ((double)rand() / (RAND_MAX));
			if (rngJesus <= this->dying) {
				Tour* eve = new Tour(toursize);
				eve->generateTour(this->map);
				newPop->saveTour(i, eve);
				continue;
			}
			Tour* a = new Tour(this->tournamentSelection(i));
			newPop->saveTour(i, a);
		}
	}

	Tour* tournamentSelection(int pos) {
		Population tourn(this->tournamentSize, this->popSize, this->map);
		tourn.saveTour(0, this->pop->getTour(pos));
		for (int z = 1; z < this->tournamentSize; z++) {
			int rPos = rand() % this->popSize;
			tourn.saveTour(z, this->pop->getTour(rPos));
		}
		Tour* fittest = tourn.getFittest();
		tourn.setNull();
		return fittest;
	}
	
	void doLogging(int gen) {
		Tour* top = this->pop->getFittest();
		if (top->getFitness() > bestTour.second->getFitness()) {
			delete this->bestTour.second;
			this->bestTour = std::make_pair(gen, new Tour(top));
			// this->stats->logTopCandidate(gen, top, this->map);
		}
		this->stats->logData(gen, pop->getAvgFitness(), this->bestTour.second);
	}

	void evolvePopulation(int gen) {
		// this->pop->checkForErrorTours();
		std::cout << "Gen: " << gen << ", Avg. Fitness: " << this->pop->calcAverageFitness() << std::endl;
		doLogging(gen);
		int toursize = this->map->size;
		Population* newPop = new Population(this->popSize, toursize,this->map);
		select(newPop,toursize);
		selectCrossOver(newPop, toursize);
		this->mutate(newPop, toursize);
		delete pop;
		this->pop = newPop;
	}
	
public:
	genAlgo() {
		pop = new Population(this->map, this->popSize, this->map->size);	
		this->stats = new Statistics();
		bestTour = std::make_pair(0,new Tour(this->pop->getTour(0)));
		// pop->printAll();
	}

	~genAlgo() {}

	void evolveTillTimesUp(int minutes) {
			std::thread t(timer,minutes);
			int curGen = 0;
			this->stats->logBasics(this->popSize, this->gens, this->mutationRate, this->crossOverRate, this->dying);
			while (!abortEvol) {
				this->evolvePopulation(curGen);
				curGen++;
			}
			t.join();
			this->stats->logTopCandidate(this->bestTour.first, this->bestTour.second, this->map);
			this->stats->logTime(minutes);
		}
};

int main(int argc, char *argv[]) {
	std::srand(unsigned(std::time(0)));
	abortEvol = false;
	//omp_set_dynamic(0);     // Explicitly disable dynamic teams
	//omp_set_num_threads(NT); // Use 4 threads for all consecutive parallel regions
	genAlgo Algo;
	// Algo.startEvolving();
	Algo.evolveTillTimesUp(15);
	return 0;
}
