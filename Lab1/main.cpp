#include "intsort.cpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex> // not necessary in 2012... may break older compilers
#include <chrono>
#include <unistd.h>
#include <fstream>

using namespace std;

//-------GLOBAL--------
typedef std::vector< std::string >  stringVector;
const size_t NUM_OF_THREADS = 8;
std::mutex globalMutex;
stringVector* lines;
vector<vector<string>*> allVecs;
int inputSize;
int tID = -1;
//-------GLOBAL--------

void threadStdSort(int iD) {
    std::cout << "hello from thread: " << iD << std::endl;
    stringVector* s = new stringVector();
    for(int i = 0; i < inputSize / NUM_OF_THREADS; i++) {
        if(iD+i*NUM_OF_THREADS >= inputSize) {
            break;
        }
        s->push_back(lines->at(iD+i*NUM_OF_THREADS));
    }
//    for(int i =0; i < s->size(); i++) {
//        std::cout << "thread "<< iD <<" gehört: " << s->at(i) << std::endl;
//    }
    std::sort(s->begin(),s->end(),[](string a, string b) {
                  return a<b;
              });
    allVecs.push_back(s);
}

void mergeSubVecs(vector<string>* v1, vector<string>* v2) {
     vector<string>* res = new vector<string>();
     int n1 = v1->size();
     int n2 = v2->size();
     int i = 0;
     int j = 0;
     while (i < n1 && j < n2) {
         if (v1->at(i) <= v2->at(j)) {
             res->push_back(v1->at(i));
             i++;
         } else {
             res->push_back(v2->at(j));
             j++;
         }
     }
     //copy remaining elements
     while (i < n1) {
         res->push_back(v1->at(i));
         i++;
     }
     while (j < n2) {
         res->push_back(v2->at(j));
         j++;
     }     
     lines = res;
     res = NULL;
     delete res;
     //allVecs.push_back(res);
 }

void normalSort() {
    std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();
    std::thread* threads[NUM_OF_THREADS];
    string inputLine;
    lines = new stringVector();
//    std::cout << "before read"<< std::endl;
//    while(!std::cin.fail()) {
//        getline( std::cin, inputLine);
//        lines->push_back(inputLine);
//    }

//   lines = new stringVector{"abc","ab","dieter","dagma","bernd","ctf","supi","acht","niemand","ende","asf","asedfasdf","sdfgs","lkhjk","ouzgz"};
//    string fname = "all_wo_country.dic";
//    string fname = "test.txt";
	string fname = "randProcessed";
    std::fstream file;
    file.open(fname,ios::in);
    if(file.is_open()) {
        std::cout << "open" << std::endl;
    } else {
        std::cout << "notOpen" << std::endl;
    }
    while(getline(file,inputLine)) {
        lines->push_back(inputLine);
    }

    std::cout << "lines size after input: "<< lines->size() << std::endl;

    inputSize = lines->size();
    for(int i=0;i<NUM_OF_THREADS;i++) {
            threads[i] = new std::thread(&threadStdSort,i);
    }
    for(int i=0;i<NUM_OF_THREADS;i++) {
        threads[i]->join();
        delete threads[i];
        threads[i] = NULL;
    }
//    for(int z=0; z < allVecs.size(); z++) {
//        std::cout << z << std::endl;
//        for(int x=0; x< allVecs[z]->size(); x++) {
//            std::cout << allVecs[z]->at(x) << std::endl;
//        }
//    }


    lines->clear();
    lines = allVecs[0];
    allVecs.erase(allVecs.begin());
    for(int z=NUM_OF_THREADS; z > 1; z--) {
        mergeSubVecs(lines, allVecs[0]);
        allVecs.erase(allVecs.begin());
    }
    printf("\n");
    ofstream res;
    res.open("res2.txt",ios::out);
    for(int i=0;i<lines->size();i++) {
        std::cout << lines->at(i) << std::endl;
    }


    std::chrono::system_clock::time_point endTime = std::chrono::system_clock::now();
    std::chrono::microseconds microRunTime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    double runTime = microRunTime.count() / 1000000.0;
	printf("\n");
    std::cout << std::setprecision( 8 )
              << "Wall clock time = " << runTime << " seconds."
              << std::endl << std::flush;
    res.close();
}

int main() {
    normalSort();
    return 0;
}

