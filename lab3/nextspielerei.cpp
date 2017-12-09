#include <iostream>
#include <string.h>
#include <complex>
#include "CMatrix.h"
#include "mpi.h"

using namespace std;

int main(int argc, char** argv) {
    MPI_Status status;
    int rank, np, peer;
    int length;
    char name[MPI_MAX_PROCESSOR_NAME + 1];
    const int msize = 4;

    int closeToRandomVariable = 1;

    int A[msize][msize] = {
        {1, 1, 2, 4},
        {3, 4, 5, 6},
        {6, 7, 8, 9},
        {0, 1, 2, 3}
    };
    int B[msize][msize] = {
        {8, 7, 6, 5},
        {5, 2, 3, 2},
        {2, 1, 4, 1},
        {0, 1, 2, 6}
    };
    int C[msize][msize];

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
	

    if (rank == 0) {
        for (int i = 0; i < msize; i++) {
            for (int j = 0; j < msize; j++) {
                int a, b;
                a = A[i][j];
                b = B[i][j];
                // cout << "I am sending: " << a << " and 
                //   << b << " to " << j + (i * msize) << endl;
                int destination = j + (i * msize);
                //int destination_j = ((j+i) + (j * msize)) % (msize * (j+1));
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

	if (rank == closeToRandomVariable) {
		cout << "MAtrix A shifted: " << endl;
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
        a = e;
        b = d;
    }

    for (int i = 0; i < msize; i++) {
        c += a*b;
        if (rank == closeToRandomVariable){
            cout << "I multiplied a*b=c " << a << "*" << b << "=" << c << endl;
		}
        int row_Dest, col_Dest;
        //(n + (i % n)) % n
        row_Dest = (row_rank - 1);
        row_Dest = (row_Dest + msize) % msize;


        col_Dest = (col_rank - 1);
        col_Dest = (col_Dest + msize) % msize;
		if (rank == closeToRandomVariable) {
			cout << "i output to ROW: " << row_Dest << " And COL: " << col_Dest <<  endl;
		}

        a2 = a;
        b2 = b;

        MPI_Sendrecv(&a, 1, MPI_INT, row_Dest, 0, &a2, 1, MPI_INT, MPI_ANY_SOURCE, 0, (row_rank + 1) % msize, &statusForA);
        MPI_Sendrecv(&b, 1, MPI_INT, col_Dest, 0, &b2, 1, MPI_INT, MPI_ANY_SOURCE, 0, (col_rank + 1) % msize, &statusForB);
        a = a2;
        b = b2;

		if (rank == closeToRandomVariable) {
			cout << "I Got A : " << a <<"from: " << (row_rank + 1) % msize << " and B:  " << b << " from: "<< (col_rank + 1) % msize << endl;
		}
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
    return 0;
}