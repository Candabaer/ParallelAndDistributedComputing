#include <iostream>
#include <string.h>
#include "CMatrix.h"
#include "mpi.h"

using namespace std;
int main(int argc, char** argv)
{
    MPI_Status  status;
    int         rank, np, peer;
    int         i, j, length;
    char name[MPI_MAX_PROCESSOR_NAME+1];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    
    
    if (MPI_SUCCESS != MPI_Get_processor_name(name, &length)) {
        strcpy(name, "unknown");
    }
    
    for (i=0; i<np; i++) {
        if (rank==i) { /* rank i talks to each higher rank */
            for(j=i+1; j<np; j++) {
                printf("checking connection between rank %d on ... %s and rank %-4d\n",i, name, j);
                MPI_Send(&rank, 1, MPI_INT, j, rank, MPI_COMM_WORLD);
                MPI_Recv(&peer, 1, MPI_INT, j, j, MPI_COMM_WORLD, &status);
            } 
        }  
        else if (rank>i) {
            /* receive from and reply to rank i */
            MPI_Recv(&peer, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Send(&rank, 1, MPI_INT, i, rank, MPI_COMM_WORLD);
        }      
    }
    
    
    
    
    
    
    if (0 == rank)
        printf("Connectivity test on %d processes PASSED.\n", np); 
    MPI_Finalize();
return 0;
 
}


