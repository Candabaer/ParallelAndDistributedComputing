#include <iostream>
#include <string.h>
#include <complex>
#include "CMatrix.h"
#include "mpi.h"
#include <vector>

using namespace std;

MPI_Status status;
int rank, np, peer;
int length;
char name[MPI_MAX_PROCESSOR_NAME + 1];
const int msize = 8;

int closeToRandomVariable = 1;

int A[msize][msize] = {
	{ 1, 1, 2, 4, 5, 6, 7, 8 },
	{ 3, 4, 5, 6, 5, 6, 7, 8 },
	{ 6, 7, 8, 9, 5, 6, 7, 8 },
	{ 0, 1, 2, 3, 4, 5, 6, 7 },
	{ 1, 1, 2, 4, 5, 6, 7, 8 },
	{ 3, 4, 5, 6, 5, 6, 7, 8 },
	{ 6, 7, 8, 9, 5, 6, 7, 8 },
	{ 0, 1, 2, 3, 4, 5, 6, 7 }
};
int B[msize][msize] = {
	{ 8, 7, 6, 5, 4, 3, 2, 1 },
	{ 5, 2, 3, 2, 1, 9, 8, 7 },
	{ 2, 1, 4, 1, 5, 5, 5, 5 },
	{ 0, 1, 2, 6, 5, 5, 5, 5 },
	{ 8, 7, 6, 5, 4, 3, 2, 1 },
	{ 5, 2, 3, 2, 1, 9, 8, 7 },
	{ 2, 1, 4, 1, 5, 5, 5, 5 },
	{ 0, 1, 2, 6, 5, 5, 5, 5 }
};
int C[msize][msize];

vector<int[][]> blocksA;
vector<int[][]> blocksB;


void writeMatrix() {

}

void makeBlocks() {
	int blockSize = msize/sqrt(np); // 2x2
	int aB = msize / blockSize;	// 4 blocks
	for (int br = 0; br < aB; br++) {
		for (int bc = 0; bc < aB; bc++) {
			int blockA[blockSize][blockSize];
			int blockB[blockSize][blockSize];
			copy(&A[br*blockSize][bc*blockSize], &A[br*blockSize][bc*blockSize] + blockSize*blockSize, &blockA[0][0]);
			copy(&B[br*blockSize][bc*blockSize], &B[br*blockSize][bc*blockSize] + blockSize*blockSize, &blockB[0][0]);
			blocksA.push_back(blockA);
			blocksB.push_back(blockB);
		}
	}
}


int main(int argc, char** argv) {
	makeBlocks();
	cout << " ABlocks: " << blocksA.size() << endl;
	cout << " BBlocks: " << blocksB.size() << endl;




	/*

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    MPI_Comm rowCom, colCom;
    MPI_Comm_split(MPI_COMM_WORLD, rank / msize, rank, &rowCom);
    MPI_Comm_split(MPI_COMM_WORLD, rank % msize, rank, &colCom);

    int row_rank, row_size, col_rank, col_size;

    MPI_Comm_rank(rowCom, &row_rank);
    MPI_Comm_size(rowCom, &row_size);

    MPI_Comm_rank(colCom, &col_rank);
    MPI_Comm_size(colCom, &col_size);

    MPI_Request requestForA, requestForB, requestForC;
    MPI_Status statusForA, statusForB, statusForC;

    int c = 0;
    int d, e;
    e = 0;
    d = 0;

    int offset = 0;
	
	int tmp;
	int TMP_A[msize][msize];
	int TMP_B[msize][msize];
	
	copy(&A[0][0], &A[0][0]+msize*msize, &TMP_A[0][0]);
	copy(&B[0][0], &B[0][0]+msize*msize, &TMP_B[0][0]);

	for (int i = 0; i<msize; i++){
		for (int j = 0; j<msize; j++){
			int n = (j-i)%msize;
			n = (msize + (n%msize)) % msize;
			int m = (i-j)%msize;
			m = (msize + (m%msize)) % msize;
			A[i][n] = TMP_A[i][j];
			B[m][j] = TMP_B[i][j];
		}
	}
	
	if (rank == closeToRandomVariable) {
		cout << "Matrix A shifted: " << endl;
		for (int i = 0; i < msize; i++) {
			for (int j = 0; j < msize; j++) {
				cout << A[i][j] << " ";
			}
			cout << endl;
		}

		cout << endl;
		cout << "Matrix B shifted: " << endl;
		for (int i = 0; i < msize; i++) {
			for (int j = 0; j < msize; j++) {
				cout << B[i][j] << " ";
			}
			cout << endl;
		}
	}

    if (rank == 0) {
        for (int i = 0; i < msize; i++) {
            for (int j = 0; j < msize; j++) {
                int a, b;
                a = A[i][j];
                b = B[i][j];
                int destination = j + (i * msize);
                if (0 == destination) {
                    d = a;
                    e = b;
                } else {
                    MPI_Send(&a, 1, MPI_INT, destination, 0, MPI_COMM_WORLD);
                    MPI_Send(&b, 1, MPI_INT, destination, 0, MPI_COMM_WORLD);
                }
            }
        }
    }

    int a, b;
    int a2, b2;
    a2 = 0;
    b2 = 0;
    if (rank != 0) {
        MPI_Recv(&a, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
       // cout << "Rank " << rank << "about to receive a whpich = " << a << endl;
        MPI_Recv(&b, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
       // cout << "Rank " << rank << "about to receive b which = " << b << endl;
    }	

    if (rank == 0) {
        a = d;
        b = e;
    }

	cerr << rank << " ,a:	" << a << ", b:	" << b << endl;

    for (int i = 0; i < msize; i++) {
        c += a*b;
        if (rank == closeToRandomVariable){
            // cout << "I multiplied a*b=c " << a << "*" << b << "=" << c <<" ,my A is: " << endl;

		}

        int row_Dest, col_Dest;
        //(n + (i % n)) % n
        row_Dest = (row_rank - 1);
        row_Dest = (row_Dest + msize) % msize;


        col_Dest = (col_rank - 1);
        col_Dest = (col_Dest + msize) % msize;

        a2 = a;
        b2 = b;

        MPI_Sendrecv(&a, 1, MPI_INT, row_Dest, 0, &a2, 1, MPI_INT, (row_rank + 1) % msize, 0, rowCom, &statusForA);
        MPI_Sendrecv(&b, 1, MPI_INT, col_Dest, 0, &b2, 1, MPI_INT, (col_rank + 1) % msize, 0, colCom, &statusForB);
        a = a2;
        b = b2;

		//if (rank == closeToRandomVariable) {
		//	cout << "I Got A: " << a << " from: " << (row_rank + 1) % msize << " and B:  " << b << " from: " << (col_rank + 1) % msize << endl;
		//}
    }

	if (rank != 0) {
		MPI_Send(&c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
    if (rank == 0) {
        for (int i = 0; i < msize; i++) {
            for (int j = 0; j < msize; j++) {
                int destination = j + (i * msize);
				if (destination != 0) {
					MPI_Recv(&c, 1, MPI_INT, destination, 0, MPI_COMM_WORLD, &status);
				}
                C[i][j] = c;
				//cout << "C[" << i << "]["<< j << "] = " << c << endl;
				cout << C[i][j] << " "; 
            }
			cout << endl;
        }
    }

    MPI_Finalize();
	*/
    return 0;
}