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

    int A[msize][msize] = {
        {0, 1, 2, 4},
        {3, 4, 5, 6},
        {6, 7, 8, 9},
        {0, 1, 2, 3}
    };
    int B[msize][msize] = {
        {8, 7, 6, 5},
        {5, 4, 3, 2},
        {2, 1, 0, 1},
        {0, 1, 2, 3}
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

    MPI_Comm_rank(rowCom, &col_rank);
    MPI_Comm_size(rowCom, &col_size);


    MPI_Request requestForA, requestForB, requestForC;
    MPI_Status statusForA, statusForB, statusForC;

    //printf("WORLD rank: %d \t  Row_RANK/Col_Rank: %d/%d\n",
    //        rank ,row_rank, col_rank);

    //printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",
    ///      rank, np, col_rank, col_size);

    int c = 0;

    //Verteilung alle werte aus beiden matrizen;
    if (rank == 0) {
        for (int i = 0; i < msize; i++) {
            for (int j = 0; j < msize; j++) {
                int a, b;
                a = A[i][j];
                b = B[i][j];
                // cout << "I am sending: " << a << " and " 
                //   << b << " to " << j + (i * msize) << endl;
                
                MPI_Isend(&a, 1, MPI_INT, j + (i * msize), 0, MPI_COMM_WORLD, &requestForA);
                //MPI_Wait(&requestForA, &statusForA);
                MPI_Isend(&b, 1, MPI_INT, j + (i * msize), 0, MPI_COMM_WORLD, &requestForB);
                if(rank == 0 && j + (i * msize) == 0)
                    cout << "habidabidubasdaGFGdasgagEG" << endl;
                if(rank == 3){
                    cout << "\t DAFAQ" << endl;
                }
                //MPI_Wait(&requestForB, &statusForB);
            }
        }
    }
    //Bis hier hin kein Deadlock
    //Speicher von Array A und B freigeben um mehr freien speicher zu haben;

    // JEDER PROZESS HAT SEINEN EIGEN VALUE IN C
    bool notEvenOnce = true;
    for (int i = 0; i < msize; i++) {
        int a, b;

        if (notEvenOnce) {
            notEvenOnce = false;
            if (rank == 0)
                cout << "FICK DICH UND HALTS MAUL" << endl;
            MPI_Recv(&a, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&b, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        } else {
            if (rank == 0)
                cout << "Maybe my life is a sad piece of shit and i should've studied\n"
                    "applied arts like my mother always wanted" << endl;
            MPI_Irecv(&a, 1, MPI_INT, MPI_ANY_SOURCE, 0, rowCom, &requestForA);
            MPI_Irecv(&b, 1, MPI_INT, MPI_ANY_SOURCE, 0, colCom, &requestForB);
            if (rank == 0)
                cout << "I've messed up badly if I reached this fucking point" << endl;
        }
        if (!notEvenOnce) {
            if (rank == 0)
                cout << "In the end I should be here: " << endl;

            MPI_Isend(&a, 1, MPI_INT, (row_rank + 1) % msize, 0, rowCom, &requestForA);
            MPI_Wait(&requestForA, &statusForA);
            MPI_Isend(&b, 1, MPI_INT, (col_rank + 1) % msize, 0, colCom, &requestForB);
            MPI_Wait(&requestForB, &statusForB);

            if (rank == 0)
                cout << "But the real question can I get here" << endl;
        }
        c += a*b;
    }

    cout << endl;
    
    cout << "Rank: " << rank << " reporting in" << endl;   
    //Hängt sich gerade auf weil 0 mit sich selber redet dieser pisser!!!!!!!!
    if (rank == 0)
        cout << "look at me I can't even not deadlock for fucking once " << endl;
    //MPI_Wait(&requestForC, &statusForC);
    /*MPI_Isend(&c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &re
     * questForC);
     */ 
    MPI_Isend(&c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &requestForC);
    MPI_Wait(&requestForC, &status);
    
    //MPI_Recv(&a, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

    if (rank == 0)
        cout << "Pls just fucking work 000000000000 " << endl;
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
                MPI_Recv(&c, 1, MPI_INT, j + (i * msize), 0, MPI_COMM_WORLD, &status);
                C[i][j] = c;
                cout << c << endl;
                cout << "faggotry should be enforced" << endl;
            }
        }
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