#include <iostream>
#include <string.h>
#include <sstream>
#include <complex>
#include "CMatrix.h"
#include "mpi.h"
#include <vector>
#include <cstdlib>

using namespace std;

MPI_Status status;
int rank, np, peer;
int length;
char name[MPI_MAX_PROCESSOR_NAME + 1];
const int totalMSize = 8;
int blockSize;
int closeToRandomVariable = 2;
int aB;
int sqrtAB;

int A[totalMSize][totalMSize] = {
	{ 65, 66, 67, 68, 69, 70, 71, 72 },
	{ 73, 10, 11, 12, 13, 12, 15, 16 },
	{ 17, 18, 19, 20, 21, 22, 23, 24 },
	{ 25, 26, 27, 28, 29, 30, 31, 32 },
	{ 33, 34, 35, 36, 37, 38, 39, 40 },
	{ 41, 42, 43, 44, 45, 46, 47, 48 },
	{ 49, 50, 51, 52, 53, 54, 55, 56 },
	{ 57, 58, 59, 60, 61, 62, 63, 64 }
};
int B[totalMSize][totalMSize] = { 
	{ 65, 66, 67, 68, 69, 70, 71, 72 },
	{ 73, 10, 11, 12, 13, 12, 15, 16 },
	{ 17, 18, 19, 20, 21, 22, 23, 24 },
	{ 25, 26, 27, 28, 29, 30, 31, 32 },
	{ 33, 34, 35, 36, 37, 38, 39, 40 },
	{ 41, 42, 43, 44, 45, 46, 47, 48 },
	{ 49, 50, 51, 52, 53, 54, 55, 56 },
	{ 57, 58, 59, 60, 61, 62, 63, 64 }
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

string to_string(int egal){
	ostringstream tmp;
	tmp << egal;
	return tmp.str();
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

//Shifts a Matrix for singleLine Values/each Process has one Element
void initialSingleLineShift() {
	int TMP_A[totalMSize][totalMSize];
	int TMP_B[totalMSize][totalMSize];

	copy(&A[0][0], &A[0][0] + totalMSize*totalMSize, &TMP_A[0][0]);
	copy(&B[0][0], &B[0][0] + totalMSize*totalMSize, &TMP_B[0][0]);
	for (int i = 0; i < totalMSize; i++) {
		for (int j = 0; j < totalMSize; j++) {
			int n = (j - i) % totalMSize;
			n = (totalMSize + (n% totalMSize)) % totalMSize;
			int m = (i - j) % totalMSize;
			m = (totalMSize + (m%totalMSize)) % totalMSize;
			A[i][n] = TMP_A[i][j];
			B[m][j] = TMP_B[i][j];
		}
	}
}

// Shifts a Matrix with respective Blocks.
void initialShift() {
	int TMP_A[totalMSize][totalMSize];
	int TMP_B[totalMSize][totalMSize];
	copy(&A[0][0], &A[0][0] + totalMSize*totalMSize, &TMP_A[0][0]);
	copy(&B[0][0], &B[0][0] + totalMSize*totalMSize, &TMP_B[0][0]);
	for (int br = 0; br < sqrtAB; br++) {
		for (int bc = 0; bc < sqrtAB; bc++) {
			int dA = br - bc;
			dA = (dA + sqrtAB) % sqrtAB;
			int dB = bc - br;
			dB = (dB + sqrtAB) % sqrtAB;
			int** bA = blockOutOfMat(TMP_A, br*blockSize, bc*blockSize, blockSize, blockSize);
			int** bB = blockOutOfMat(TMP_B, br*blockSize, bc*blockSize, blockSize, blockSize);
			blockIntoMat(bA, A, br*blockSize, dA*blockSize, blockSize, blockSize);
			blockIntoMat(bB, B, dB*blockSize, bc*blockSize, blockSize, blockSize);
		}
	}
}

int main(int argc, char** argv) {
//-----------------------Init Part------------------------//
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	blockSize = totalMSize / sqrt(np);
	aB = (totalMSize*totalMSize) / (blockSize*blockSize);
	sqrtAB = sqrt(aB);
	MPI_Comm rowCom, colCom;
	MPI_Comm_split(MPI_COMM_WORLD, rank / (int)sqrt(aB), rank, &rowCom);
	MPI_Comm_split(MPI_COMM_WORLD, rank % (int)sqrt(aB), rank, &colCom);
	int row_rank, row_size, col_rank, col_size;
	MPI_Comm_rank(rowCom, &row_rank);
	MPI_Comm_size(rowCom, &row_size);
	MPI_Comm_rank(colCom, &col_rank);
	MPI_Comm_size(colCom, &col_size);

	MPI_Request requestForA, requestForB, requestForC;
	MPI_Status statusForA, statusForB, statusForC;

	if (rank == 0) {
		cout << "NP: " << np << endl;
		cout << "BlockSize: " << blockSize << endl;
		cout << "amountBlocks: " << aB << endl;
		initialShift();

		cout << "AShifted: " << endl;
		for (int r = 0; r < totalMSize; r++) {
			for (int c = 0; c < totalMSize; c++) {
				cout << A[r][c] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "BShifted: " << endl;
		for (int r = 0; r < totalMSize; r++) {
			for (int c = 0; c < totalMSize; c++) {
				cout << B[r][c] << " ";
			}
			cout << endl;
		}
	}

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
		int** saveA = alloc_2d_int(blockSize, blockSize);
		int** saveB = alloc_2d_int(blockSize, blockSize);
		for (int br = 0; br < sqrtAB; br++) {
			for (int bc = 0; bc < sqrtAB; bc++) {
				int destination = bc + (br * sqrt(aB));
				BlockA = blockOutOfMat(A, br*blockSize, bc*blockSize, blockSize, blockSize);
				BlockB = blockOutOfMat(B, br*blockSize, bc*blockSize, blockSize, blockSize);
				printMatrix(BlockB, blockSize, blockSize, to_string(destination));
				if (0 == destination) {
					std::copy(&BlockA[0][0],&BlockA[0][0]+blockSize*blockSize,&saveA[0][0]);
					std::copy(&BlockB[0][0],&BlockB[0][0]+blockSize*blockSize,&saveB[0][0]);
				}
				else {
					MPI_Send(&BlockA[0][0], blockSize*blockSize, MPI_INT, destination, 0, MPI_COMM_WORLD);
					MPI_Send(&BlockB[0][0], blockSize*blockSize, MPI_INT, destination, 0, MPI_COMM_WORLD);
					// cout << "Test after sending blocks!" << endl;
				}
			}
		}
		std::copy(&saveA[0][0],&saveA[0][0]+blockSize*blockSize,&BlockA[0][0]);
		std::copy(&saveB[0][0],&saveB[0][0]+blockSize*blockSize,&BlockB[0][0]);
	//	freeMem(saveA, blockSize);
	//	freeMem(saveB, blockSize);
	}
//-----------------------DO MPI WITH BLOCKS------------------------//
	int** a2 = alloc_2d_int(blockSize, blockSize);
	int** b2 = alloc_2d_int(blockSize, blockSize);
	if (rank != 0) {
		MPI_Recv(&BlockA[0][0], blockSize*blockSize, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		// cout << "Rank " << rank << "about to receive a whpich = " << a << endl;
		MPI_Recv(&BlockB[0][0], blockSize* blockSize, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		// cout << "Rank " << rank << "about to receive b which = " << b << endl;
	}
	// cerr << rank << " ,a:	" << a << ", b:	" << b << endl;
	for (int i = 0; i < sqrtAB; i++) {
		if (rank == closeToRandomVariable) {
			printMatrix(BlockC, blockSize, blockSize, "Block C before multAdd");
		}
		addMult(BlockA, BlockB, blockSize, BlockC);
		if (rank == closeToRandomVariable) {
			cout << "I multiplied a*b=c " << endl;
			printMatrix(BlockA, blockSize, blockSize, "BlockA");
			printMatrix(BlockB, blockSize, blockSize, "BlockB");
			printMatrix(BlockC, blockSize, blockSize, "BlockC");
		}
		int row_Dest, col_Dest;
		row_Dest = (row_rank - 1);
		row_Dest = (row_Dest + sqrtAB) % sqrtAB;
		col_Dest = (col_rank - 1);
		col_Dest = (col_Dest + sqrtAB) % sqrtAB;
		//std::copy(&BlockA[0][0],&BlockA[0][0]+blockSize*blockSize,&a2[0][0]);
		//std::copy(&BlockB[0][0],&BlockB[0][0]+blockSize*blockSize,&b2[0][0]);
		if (rank == closeToRandomVariable) {
			cout << "RowDest: " << row_Dest << " , row from: " << (row_rank + 1) % sqrtAB << endl;
			cout << "colDest: " << col_Dest << " , row from: " << (col_rank + 1) % sqrtAB << endl;
		}
		MPI_Sendrecv(&BlockA[0][0], blockSize*blockSize, MPI_INT, row_Dest, 0, &a2[0][0], blockSize*blockSize, MPI_INT, (row_rank + 1) % sqrtAB, 0, rowCom, &statusForA);
		MPI_Sendrecv(&BlockB[0][0], blockSize*blockSize, MPI_INT, col_Dest, 0, &b2[0][0], blockSize*blockSize, MPI_INT, (col_rank + 1) % sqrtAB, 0, colCom, &statusForB);
		std::copy(&a2[0][0],&a2[0][0]+blockSize*blockSize,&BlockA[0][0]);
		std::copy(&b2[0][0],&b2[0][0]+blockSize*blockSize,&BlockB[0][0]);
	}
	if (rank != 0) {
		MPI_Send(&BlockC[0][0], blockSize*blockSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
//----------------------------Recv Blocks------------------------//
	if (rank == 0) {
		for (int i = 0; i < sqrtAB; i++) {
			for (int j = 0; j < sqrtAB; j++) {
				int destination = j + (i * sqrtAB);
				if (destination != 0) {
					MPI_Recv(&BlockC[0][0], blockSize*blockSize, MPI_INT, destination, 0, MPI_COMM_WORLD, &status);
				}
				blockIntoMat(BlockC, C, i*blockSize, j*blockSize, blockSize, blockSize);
			}
		}
		for (int r = 0; r < totalMSize; r++) {
			for (int c = 0; c < totalMSize; c++) {
				cout << C[r][c] << " ";
			}
			cout << endl;
		}	
	}
	MPI_Finalize();	
    return 0;
}
