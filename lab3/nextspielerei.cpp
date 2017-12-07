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

    int closeToRandomVariable = 6;

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

    /*Communicator definiert*/

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

    //printf("WORLD rank: %d \t  Row_RANK/Col_Rank: %d/%d\n",
    //        rank ,row_rank, col_rank);

    //printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",
    ///      rank, np, col_rank, col_size);

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

    //Verteilung alle werte aus beiden matrizen;
    if (rank == 0) {
        for (int i = 0; i < msize; i++) {
            for (int j = 0; j < msize; j++) {
                int a, b;
                a = A[i][j];
                b = B[i][j];
                // cout << "I am sending: " << a << " and 
                //   << b << " to " << j + (i * msize) << endl

                int destination = j + (i * msize);
                //int destination_j = ((j+i) + (j * msize)) % (msize * (j+1));
                if (0 == destination) {
                    d = a;
                    e = b;
                } else {
                    MPI_Send(&a, 1, MPI_INT, destination, 0, MPI_COMM_WORLD);
                    MPI_Send(&b, 1, MPI_INT, destination, 0, MPI_COMM_WORLD);
                    //cout << "Outputting a & b: " << a << "     " <<  b << endl;
                }

                //cout << "Rank: " << rank << " send a: " << a << " with b: " 
                //        << b << " to " << destination << endl;

                //MPI_Wait(&requestForB, &statusForB);
            }
        }
    }
    //Bis hier hin kein Deadlock
    //Speicher von Array A und B freigeben um mehr freien speicher zu haben;

    // JEDER PROZESS HAT SEINEN EIGEN VALUE IN C

    /**Entwicklertagebuch 
     ich muss bei der ersten iteration ein bsishen was verschieben
     * 
     * Ey kinder was macht ihr da aufm dach struggeln & husseln
     mussdeshalb leider 1mal eine unnötige kommunikation machen.
     */

    int a, b;
    int a2, b2;
    a2 = 0;
    b2 = 0;
    if (rank != 0) {
        MPI_Recv(&a, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        //cout << "Rank " << rank << "about to receive a whpich = " << a << endl;
        MPI_Recv(&b, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        //cout << "Rank " << rank << "about to receive b which = " << b << endl;
    }
    if (rank == 0) {
        a = e;
        b = d;
    }


    //First time shifting damit alles seine richtigkeit hat
    bool beGentleSenpai = true;

    /*if (beGentleSenpai) {
        int row_Dest, col_Dest;
        //(n + (i % n)) % n
        row_Dest = (row_rank) % msize;
        row_Dest = (msize + (row_Dest % msize)) % msize;

        col_Dest = (col_rank) % msize;
        col_Dest = (msize + (col_Dest % msize)) % msize;

        a2 = a;
        b2 = b;
        int offset_i = 0;
        int offset_j = 0;
        if(rank % msize){
            offset_i = rank / msize;
        }
        MPI_Sendrecv(&offset_i, 1, MPI_INT, row_Dest, 0, &offset_i_, 1, MPI_INT, MPI_ANY_SOURCE, 0, rowCom, &statusForA);
        MPI_Sendrecv(&offset_j, 1, MPI_INT, row_Dest, 0, &offset_j_, 1, MPI_INT, MPI_ANY_SOURCE, 0, rowCom, &statusForA);
        
        //if(rank / msize == )
        row_Dest %= (row_Dest + offset_i_) % msize;
        col_Dest %= (col_Dest + offset_j_) % msize;
        MPI_Sendrecv(&a, 1, MPI_INT, row_Dest, 0, &a2, 1, MPI_INT, MPI_ANY_SOURCE, 0, rowCom, &statusForA);
        MPI_Sendrecv(&b, 1, MPI_INT, col_Dest, 0, &b2, 1, MPI_INT, MPI_ANY_SOURCE, 0, colCom, &statusForB);

        a = a2;
        b = b2;
        beGentleSenpai = false;
    }*/


    for (int i = 0; i < msize; i++) {

        // if(rank == 0 )
        //    cout << "E = " << e << " D = " << d << endl; 

        c += a*b;
        if (rank == closeToRandomVariable){
            //cout << "I multiplied a*b=c " << a << "*" << b << "=" << c << endl;
	}
        int row_Dest, col_Dest;
        //(n + (i % n)) % n
        row_Dest = (row_rank - 1) % msize;
        row_Dest = (msize + (row_Dest % msize)) % msize;

        col_Dest = (col_rank - 1) % msize;
        col_Dest = (msize + (col_Dest % msize)) % msize;

        //cout << "Row_Rank: " << row_Dest << endl;


        a2 = a;
        b2 = b;
        //cout <<  " Rank "<< rank <<  " Row_Rank: " << row_rank << " Col_Rank: " << col_rank << " sends: "
        //        << "\t A = " << a << " B = " << b << endl;
        MPI_Sendrecv(&a, 1, MPI_INT, row_Dest, 0, &a2, 1, MPI_INT, MPI_ANY_SOURCE, 0, rowCom, &statusForA);
        MPI_Sendrecv(&b, 1, MPI_INT, col_Dest, 0, &b2, 1, MPI_INT, MPI_ANY_SOURCE, 0, colCom, &statusForB);
        a = a2;
        b = b2;

        if (rank == closeToRandomVariable)
            cout << "Outputting a : " << a << "   B:  " << b << endl;


        /*MPI_Send(&a, 1, MPI_INT, row_Dest, 0, rowCom);
        //cout << "ROW_Rank " << row_rank << "about to send a to" << row_Dest  << endl;
        MPI_Send(&b, 1, MPI_INT, col_Dest, 0, colCom);
        //cout << "COL_Rank " << col_rank << "about to send b to" << col_Dest  << endl;
        

        MPI_Recv(&a, 1, MPI_INT, MPI_ANY_SOURCE, 0, rowCom, &statusForA);
        cout << "ROW " <<row_rank << "about to rec a" << endl;
        MPI_Recv(&b, 1, MPI_INT, MPI_ANY_SOURCE, 0, colCom, &statusForB);
        cout << "COL " <<col_rank << "about to rec b" << endl;
         */
    }


    //   cout << "Rank: " << rank << " reporting in" << endl;   
    //Hängt sich gerade auf weil 0 mit sich selber redet dieser pisser!!!!!!!!
    //if (rank == 0)
    //cout << "look at me I can't even not deadlock for once " << endl;
    //MPI_Wait(&requestForC, &statusForC);
    /*MPI_Isend(&c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &re
     * questForC);
     */
    if (rank != 0)
        MPI_Send(&c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    //MPI_Recv(&a, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

    //if (rank == 1)
    //    cout << "cmon bruh 1111111111" << endl;  
    //if (rank == 0)
    //  cout << "Hier komm ich aber schon an" << endl;
    if (rank == 0) {
        for (int i = 0; i < msize; i++) {
            for (int j = 0; j < msize; j++) {
                //MPI_Wait(&requestForC, &status);
                //MPI_Irecv(&c, 1, MPI_INT, j + (i * msize), 0, MPI_COMM_WORLD,
                //        &requestForC);
                //MPI_Wait(&requestForC, &status);
                int destination = j + (i * msize);
                if (destination != 0)
                    MPI_Recv(&c, 1, MPI_INT, destination, 0, MPI_COMM_WORLD, &status);

                C[i][j] = c;

         //cout << "C[" << i << "]["<< j << "] = " << c << endl;
		cout << C[i][j] << " "; 
            }
	cout << endl;
        }

        //cout << "c: " << C[2][1] << endl;
    }
    /* if (rank == 0)
         cout << "kill me pls daddy " << endl;
     if (rank == 1)
         cout << "kill me pls daddy " << endl;
     */
    /*if (rank == 0) {
        for (int i = 0; i < msize; i++) {
            cout << endl;
            for (int j = 0; j < msize; j++) {
                cout << C[i][j] << "  ";
            }

        }
    }*/

    MPI_Finalize();
    return 0;
}



/**
    
if (MPI_SUCCESS != MPI_Get_processor_name(name, &length)) {
    strcpy(name, "unknown");
}
int j=0;
int y = 3;
    
for (int i = 0; i < np; i++){
    for (int j = 0; j < np; j++){
        
        if(rank == 0){
            MPI_Send(&A[i], 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
            MPI_Send(&y, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
        }


        else if(rank == 1){
            MPI_Recv(&j, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
            MPI_Recv(&y, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
            int x = j;
            int z = y;

            cout  << "jadajada: " << x+z << endl;
        }
    }
}
 */
