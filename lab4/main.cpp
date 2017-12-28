#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


int NP = 8;

using namespace std;

int strToInt(const string &str) {
	std::stringstream temp(str);
	int res = 0;
	temp >> res;
	return res;
}

class City {
public:
	City(int ID, string name, string plz, double longi, double lat) {
		this->ID = ID;
		this->name = name;
		this->plz = plz;
		this->longi = longi;
		this->lat = lat;
		this->visited = false;
	}

private:
	int ID;
	string name;
	string plz;
	double longi;
	double lat;
	bool visited;
};

class Map {
public:
	Map(string fname) {
		ifstream is(fname);
		if (!is) {
			cerr << "Error, could not open file: " << fname << endl;
		}
		int id = 0;
		double longi, lat = 0.;
		string name, plz = "";
		while (is >> id >> plz >> longi >> lat >> name) {
			City* city = new City(id, name, plz, longi, lat);
			nodes.push_back(city);
		}
		is.close();
		cout << "Added: " << nodes.size() << " Citys." << endl;
	}

	~Map() {
	}
private:
	vector<City*> nodes;

};

int main(int argc, char *argv[]) {
	//if (argc > 3) {
	//	cerr << "Too Many Param" << endl;
	//	return -1;
	//}
	//NP = strToInt(argv[2]);
	Map* map = new Map("test_input.txt");	// (argv[1]);

}