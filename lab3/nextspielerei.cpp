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
const int totalMSize = 8;
const int blockSize = totalMSize / sqrt(np);
int closeToRandomVariable = 1;
int aB = (totalMSize*totalMSize) / (blockSize*blockSize);

int A[totalMSize][totalMSize] = {
	{ 1, 1, 2, 4, 5, 6, 7, 8 },
	{ 3, 4, 5, 6, 5, 6, 7, 8 },
	{ 6, 7, 8, 9, 5, 6, 7, 8 },
	{ 0, 1, 2, 3, 4, 5, 6, 7 },
	{ 1, 1, 2, 4, 5, 6, 7, 8 },
	{ 3, 4, 5, 6, 5, 6, 7, 8 },
	{ 6, 7, 8, 9, 5, 6, 7, 8 },
	{ 0, 1, 2, 3, 4, 5, 6, 7 }
};
int B[totalMSize][totalMSize] = {
	{ 8, 7, 6, 5, 4, 3, 2, 1 },
	{ 5, 2, 3, 2, 1, 9, 8, 7 },
	{ 2, 1, 4, 1, 5, 5, 5, 5 },
	{ 0, 1, 2, 6, 5, 5, 5, 5 },
	{ 8, 7, 6, 5, 4, 3, 2, 1 },
	{ 5, 2, 3, 2, 1, 9, 8, 7 },
	{ 2, 1, 4, 1, 5, 5, 5, 5 },
	{ 0, 1, 2, 6, 5, 5, 5, 5 }
};
int C[totalMSize][totalMSize];

int** BlockA;
int** BlockB;
int** BlockC;

void printMatrix(int** matrix, int row, int col, const string& msg) {
	std::cout << msg << std::endl;
	for (int r = 0; r < row; r++) {
		for (int c = 0; c < col; c++) {
			cout << matrix[r][c] << " ";
		}
		cout << endl;
	}
}

//Use to allocate continously memory.
int **alloc_2d_int(int rows, int cols) {
	int *data = (int *)malloc(rows*cols * sizeof(int));
	int **array = (int **)malloc(rows * sizeof(int*));
	for (int i = 0; i<rows; i++)
		array[i] = &(data[cols*i]);

	return array;
}

void freeMem(int** arr, int size) {
	for (int i = 0; i < size; ++i) {
		delete[] arr[i];//deletes an inner array of integer;
	}
	delete[] arr;
}

void initialShift() {
	int TMP_A[blockMSize][blockMSize];
	int TMP_B[blockMSize][blockMSize];
	copy(&A[0][0], &A[0][0] + blockMSize*blockMSize, &TMP_A[0][0]);
	copy(&B[0][0], &B[0][0] + blockMSize*blockMSize, &TMP_B[0][0]);
	for (int i = 0; i<blockMSize; i++) {
		for (int j = 0; j<blockMSize; j++) {
			int n = (j - i) % blockMSize;
			n = (blockMSize + (n%blockMSize)) % blockMSize;
			int m = (i - j) % blockMSize;
			m = (blockMSize + (m%blockMSize)) % blockMSize;
			A[i][n] = TMP_A[i][j];
			B[m][j] = TMP_B[i][j];
		}
	}
}

void addMult(int** mA, int** mB, int size, int** res) {
	for (int rRes = 0; rRes < size; rRes++) {
		for (int r = 0; r < size; r++) {
			for (int c = 0; c < size; c++) {
				res[rRes][r] += mA[rRes][c] * mB[c][r];
			}
		}
	}
}

int** blockOutOfMat(int m[][totalMSize], int sRow, int sCol, int dRow, int dCol) {
	int** ppBlock = alloc_2d_int(blockSize, blockSize);
	for (int r = 0, y = sRow; y < sRow + dRow; y++, r++) {
		for (int c = 0, x = sCol; x < sCol + dCol; x++, c++) {
			ppBlock[r][c] = m[y][x];
		}
	}
	return ppBlock;
}

void blockIntoMat(int** block, int res[][totalMSize], int sRow, int sCol, int dRow, int dCol) {
	for (int r = 0, y = sRow; y < sRow + dRow; y++, r++) {
		for (int c = 0, x = sCol; x < sCol + dCol; x++, c++) {
			res[y][x] = block[r][c];
		}
	}
}

int main(int argc, char** argv) {
//-----------------------Init Part------------------------//
	initialShift();
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm rowCom, colCom;
	MPI_Comm_split(MPI_COMM_WORLD, rank / aB, rank, &rowCom);
	MPI_Comm_split(MPI_COMM_WORLD, rank % aB, rank, &colCom);
	//	int row_rank, row_size, col_rank, col_size;
	MPI_Comm_rank(rowCom, &row_rank);
	MPI_Comm_size(rowCom, &row_size);
	MPI_Comm_rank(colCom, &col_rank);
	MPI_Comm_size(colCom, &col_size);

	MPI_Request requestForA, requestForB, requestForC;
	MPI_Status statusForA, statusForB, statusForC;

	BlockA = alloc_2d_int(blockSize, blockSize);
	BlockB = alloc_2d_int(blockSize, blockSize);
	BlockC = alloc_2d_int(blockSize, blockSize);
	// init all C's with 0.
	for (int r = 0; r < blockSize; r++) {
		for (int c = 0; c < blockSize; c++) {
			BlockC[r][c] = 0;
		}
	}
//-------------------Send Blocks------------------------//
	if (rank == 0) {
		int** saveA = NULL; // = alloc_2d_int(blockSize, blockSize);
		int** saveB = NULL; // = alloc_2d_int(blockSize, blockSize);
		for (int br = 0; br < sqrt(aB); br++) {
			for (int bc = 0; bc < sqrt(aB); bc++) {
				int destination = bc + (br * aB);
				BlockA = blockOutOfMat(A, br*blockSize, bc*blockSize, blockSize, blockSize);
				BlockB = blockOutOfMat(B, br*blockSize, bc*blockSize, blockSize, blockSize);
				// printMatrix(BlockA, blockSize, blockSize, "blocks");
				if (0 == destination) {
					saveA = BlockA;
					saveB = BlockB;
				}
				else {
					MPI_Send(&BlockA[0][0], blockSize*blockSize, MPI_INT, destination, 0, MPI_COMM_WORLD);
					MPI_Send(&BlockB[0][0], blockSize*blockSize, MPI_INT, destination, 0, MPI_COMM_WORLD);
					cout << "Test after sending blocks!" << endl;
				}
			}
		}
		BlockA = saveA;
		BlockB = saveB;
	}
//-----------------------DO MPI WITH BLOCKS------------------------//
	int** a2, int** b2;
	a2 = alloc_2d_int(blockSize, blockSize);
	b2 = alloc_2d_int(blockSize, blockSize);
	if (rank != 0) {
		MPI_Recv(&BlockA[0][0], blockSize*blockSize, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		// cout << "Rank " << rank << "about to receive a whpich = " << a << endl;
		MPI_Recv(&BlockB[0][0], blockSize* blockSize, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		// cout << "Rank " << rank << "about to receive b which = " << b << endl;
	}
	// cerr << rank << " ,a:	" << a << ", b:	" << b << endl;
	for (int i = 0; i < aB; i++) {
		addMult(BlockA, BlockB, blockSize, BlockC);
		if (rank == closeToRandomVariable) {
			// cout << "I multiplied a*b=c " << a << "*" << b << "=" << c <<" ,my A is: " << endl;
		}
		int row_Dest, col_Dest;
		//(n + (i % n)) % n
		row_Dest = (row_rank - 1);
		row_Dest = (row_Dest + blockSize) % blockSize;
		col_Dest = (col_rank - 1);
		col_Dest = (col_Dest + blockSize) % blockSize;
		a2 = a;
		b2 = b;
		MPI_Sendrecv(&BlockA[0][0], blockSize*blockSize, MPI_INT, row_Dest, 0, &a2[0][0], blockSize*blockSize, MPI_INT, (row_rank + 1) % aB, 0, rowCom, &statusForA);
		MPI_Sendrecv(&BlockB[0][0], blockSize*blockSize, MPI_INT, col_Dest, 0, &b2[0][0], blockSize*blockSize, MPI_INT, (col_rank + 1) % aB, 0, colCom, &statusForB);
		a = a2;
		b = b2;
	}
	if (rank != 0) {
		MPI_Send(&BlockC[0][0], blockSize*blockSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
//----------------------------Recv Blocks------------------------//
	if (rank == 0) {
		for (int i = 0; i < sqrt(aB); i++) {
			for (int j = 0; j < sqrt(aB); j++) {
				int destination = j + (i * aB);
				if (destination != 0) {
					MPI_Recv(&BlockC[0][0], blockSize*blockSize, MPI_INT, destination, 0, MPI_COMM_WORLD, &status);
				}
				blockIntoMat(BlockC, C, br*blockSize, bc*blockSize, blockSize, blockSize);
			}
		}
		for (int r = 0; r < row; r++) {
			for (int c = 0; c < col; c++) {
				cout << C[r][c] << " ";
			}
			cout << endl;
		}	
	}
	MPI_Finalize();	
    return 0;
}