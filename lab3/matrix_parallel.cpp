// file...: matrix_parallel.cpp
// desc...: parallel matrix multiplication solution -- TEMPLATE --
// <type_date_here> | <type_author_here>

#include <iostream>
#include <assert.h>
#include <mpi.h>
#include "CMatrix.h"

using namespace std;

void printError(const char* progname, const char* error) {
    if(error != NULL) {
        cerr << "ERROR: " << error << endl;
    }
    cerr << "usage: " << progname << " <matrix1> <matrix2>" << endl
    << "\twhere <matrix1> and <matrix2> are file names containing matrices." << endl;
}



// +++ main starts here +++
int main(int argc, char** argv) {

    // process arguments
    if(argc != 3) {
        printError(argv[0], "wrong number of arguments.");
        return -1;
    }
    int size, rank;
    
    
    CMatrix m1(argv[1]);                    // read 1st matrix
    CMatrix m2(argv[2]);                    // read 2nd matrix

    assert(m1.width == m2.height);          // check compatibility
    assert(m2.width == m1.height);          // check compatibility
    CMatrix result(m1.height,m2.width);     //ALLOCATE MEMORY   
    
    
    MPI::Init();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double einDimensionaleMatrixAundB[m1.width*m2.height];

    // TODO multiply matrices
    cout << "someday i will multiply here; i'm rank " << MPI::COMM_WORLD.Get_rank()
        << " of " << MPI::COMM_WORLD.Get_size() << endl;
    
    //TODO Zeilen von A einzeln aufdrößeln
    if(rank == 0){
        
    }
    
    if(rank != 0){
        
        for (int i = 0; i < m1.width; i++){
            for (int j = 0; j < m2.height; j++){
                
            }
        }
        
    }
          
    
    //TODO Spalten von B einzeln audrößeln
    
    //TODO Jedem Prozessor eine Spalte und eine Zeile geben
    
    //TODO Alles ausrechnen und in Matrix C ballern
    

    MPI::Finalize();

    // print matrix
    result.print();

    return 0;
}

// EOF
