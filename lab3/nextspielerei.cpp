#include <chrono>
#include <iostream>
#include <string.h>
#include <sstream>
#include <complex>
#include "mpi.h"
#include <vector>
#include <cstdlib>

using namespace std;

MPI_Status status;
int pID, np, peer;
int length;
char name[MPI_MAX_PROCESSOR_NAME + 1];

//int A[totalMSize][totalMSize] = {
//	{ 65, 66, 67, 68, 69, 70, 71, 72 },
//	{ 73, 10, 11, 12, 13, 14, 15, 16 },
//	{ 17, 18, 19, 20, 21, 22, 23, 24 },
//	{ 25, 26, 27, 28, 29, 30, 31, 32 },
//	{ 33, 34, 35, 36, 37, 38, 39, 40 },
//	{ 41, 42, 43, 44, 45, 46, 47, 48 },
//	{ 49, 50, 51, 52, 53, 54, 55, 56 },
//	{ 57, 58, 59, 60, 61, 62, 63, 64 }
//};
//int B[totalMSize][totalMSize] = { 
//	{ 65, 66, 67, 68, 69, 70, 71, 72 },
//	{ 73, 10, 11, 12, 13, 14, 15, 16 },
//	{ 17, 18, 19, 20, 21, 22, 23, 24 },
//	{ 25, 26, 27, 28, 29, 30, 31, 32 },
//	{ 33, 34, 35, 36, 37, 38, 39, 40 },
//	{ 41, 42, 43, 44, 45, 46, 47, 48 },
//	{ 49, 50, 51, 52, 53, 54, 55, 56 },
//	{ 57, 58, 59, 60, 61, 62, 63, 64 }
//};
//int C[totalMSize][totalMSize];

double** BlockA;
double** BlockB;
double** BlockC;

void printMatrix(double** matrix, int row, int col, const string& msg) {
	std::cout << msg << std::endl;
	for (int r = 0; r < row; r++) {
		for (int c = 0; c < col; c++) {
			cout << matrix[r][c] << " ";
		}
		cout << endl;
	}
}

int strToInt(const string &str) {
	std::stringstream temp(str);
	int res = 0;
	temp >> res;
	return res;
}
//Use to allocate continously memory.
double **alloc_2d_int(int rows, int cols) {
	double *data = (double *)malloc(rows*cols * sizeof(double));
	double **array = (double **)malloc(rows * sizeof(double*));
	for (int i = 0; i<rows; i++)
		array[i] = &(data[cols*i]);

	return array;
}

void freeMem(double** arr) {
	free(arr[0]);
	free(arr);
}

void addMult(double** mA, double** mB, int size, double** res) {
	for (int rRes = 0; rRes < size; rRes++) {
		for (int r = 0; r < size; r++) {
			for (int c = 0; c < size; c++) {
				res[rRes][r] += mA[rRes][c] * mB[c][r];
			}
		}
	}
}

double** blockOutOfMat(double** m, int sRow, int sCol, int dRow, int dCol, int blockSize) {
	double** ppBlock = alloc_2d_int(blockSize, blockSize);
	for (int r = 0, y = sRow; y < sRow + dRow; y++, r++) {
		for (int c = 0, x = sCol; x < sCol + dCol; x++, c++) {
			ppBlock[r][c] = m[y][x];
		}
	}
	return ppBlock;
}

void blockIntoMat(double** block, double** res, int sRow, int sCol, int dRow, int dCol) {
	for (int r = 0, y = sRow; y < sRow + dRow; y++, r++) {
		for (int c = 0, x = sCol; x < sCol + dCol; x++, c++) {
			res[y][x] = block[r][c];
		}
	}
}

void createRandomDouble(double& number) {
	number = rand();
	number = number / rand() * 30.14;
	if (rand() % 2) {
		number *= -1;
	}
}

int main(int argc, char** argv) {
//-----------------------Init Part------------------------//
	srand((unsigned)time(0));
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pID);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	int totalMSize;
	int closeToRandomVariable = 2;
	int blockSize;
	int aB;
	int sqrtAB;
	double** A;
	double** B;
	std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();
	if (pID == 0) {
		if (argc != 2) {
			cerr << "Not enough Arguments should be : int" << endl;
			return -1;
		}
		int mSize = strToInt(argv[1]);
		A = alloc_2d_int(mSize,mSize);
		B = alloc_2d_int(mSize, mSize);
		for (int r = 0; r < mSize; r++) {
			for (int c = 0; c < mSize; c++) {
				createRandomDouble(A[r][c]);
				createRandomDouble(B[r][c]);
			}
		}
		totalMSize = mSize;
		printMatrix(A, totalMSize, totalMSize, "Matrix A:");
		printMatrix(B, totalMSize, totalMSize, "Matrix B:");
		startTime = std::chrono::system_clock::now();
		for (int z = 1; z < np; z++) {
			MPI_Send(&mSize, 1, MPI_INT, z, 0, MPI_COMM_WORLD);
		}
	}
	else {
		MPI_Recv(&totalMSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	}
	blockSize = totalMSize / sqrt(np);
	aB = (totalMSize*totalMSize) / (blockSize*blockSize);
	sqrtAB = sqrt(aB);
	MPI_Comm rowCom, colCom;
	MPI_Comm_split(MPI_COMM_WORLD, pID / (int)sqrt(aB), pID, &rowCom);
	MPI_Comm_split(MPI_COMM_WORLD, pID % (int)sqrt(aB), pID, &colCom);
	int row_rank, row_size, col_rank, col_size;
	MPI_Comm_rank(rowCom, &row_rank);
	MPI_Comm_size(rowCom, &row_size);
	MPI_Comm_rank(colCom, &col_rank);
	MPI_Comm_size(colCom, &col_size);

	MPI_Request requestForA, requestForB, requestForC;
	MPI_Status statusForA, statusForB, statusForC;
	// Shifts a Matrix with respective Blocks.
	BlockA = alloc_2d_int(blockSize, blockSize);
	BlockB = alloc_2d_int(blockSize, blockSize);
	BlockC = alloc_2d_int(blockSize, blockSize);
	// init all C's with 0 because of += multiplication.
	for (int r = 0; r < blockSize; r++) {
		for (int c = 0; c < blockSize; c++) {
			BlockC[r][c] = 0;
		}
	}
	if (pID == 0) {
	/*	cout << "NP: " << np << endl;
		cout << "BlockSize: " << blockSize << endl;
		cout << "amountBlocks: " << aB << endl;	*/	
		double** TMP_A = alloc_2d_int(totalMSize, totalMSize);
		double** TMP_B = alloc_2d_int(totalMSize, totalMSize);
		copy(&A[0][0], &A[0][0] + totalMSize*totalMSize, &TMP_A[0][0]);
		copy(&B[0][0], &B[0][0] + totalMSize*totalMSize, &TMP_B[0][0]);
		for (int br = 0; br < sqrtAB; br++) {
			for (int bc = 0; bc < sqrtAB; bc++) {
				int dA = bc - br;
				dA = (dA + sqrtAB) % sqrtAB;
				int dB = br - bc;
				dB = (dB + sqrtAB) % sqrtAB;
				double** bA = blockOutOfMat(TMP_A, br*blockSize, bc*blockSize, blockSize, blockSize,blockSize);
				double** bB = blockOutOfMat(TMP_B, br*blockSize, bc*blockSize, blockSize, blockSize,blockSize);
				blockIntoMat(bA, A, br*blockSize, dA*blockSize, blockSize, blockSize);
				blockIntoMat(bB, B, dB*blockSize, bc*blockSize, blockSize, blockSize);
			}
		}
		freeMem(TMP_A);
		freeMem(TMP_B);
		//-------------------Send Blocks------------------------//
		double** saveA = alloc_2d_int(blockSize, blockSize);
		double** saveB = alloc_2d_int(blockSize, blockSize);
		for (int br = 0; br < sqrtAB; br++) {
			for (int bc = 0; bc < sqrtAB; bc++) {
				int destination = bc + (br * sqrt(aB));
				BlockA = blockOutOfMat(A, br*blockSize, bc*blockSize, blockSize, blockSize,blockSize);
				BlockB = blockOutOfMat(B, br*blockSize, bc*blockSize, blockSize, blockSize,blockSize);
				//printMatrix(BlockB, blockSize, blockSize, to_string(destination));
				if (0 == destination) {
					std::copy(&BlockA[0][0],&BlockA[0][0]+blockSize*blockSize,&saveA[0][0]);
					std::copy(&BlockB[0][0],&BlockB[0][0]+blockSize*blockSize,&saveB[0][0]);
				}
				else {
					MPI_Send(&BlockA[0][0], blockSize*blockSize, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD);
					MPI_Send(&BlockB[0][0], blockSize*blockSize, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD);
				}
			}
		}
		std::copy(&saveA[0][0],&saveA[0][0]+blockSize*blockSize,&BlockA[0][0]);
		std::copy(&saveB[0][0],&saveB[0][0]+blockSize*blockSize,&BlockB[0][0]);
		freeMem(saveA);
		freeMem(saveB);
		//cout << "Test after sending blocks!" << endl;
	}
//-----------------------DO MPI WITH BLOCKS------------------------//
	{
		double** a2 = alloc_2d_int(blockSize, blockSize);
		double** b2 = alloc_2d_int(blockSize, blockSize);
		if (pID != 0) {
			MPI_Recv(&BlockA[0][0], blockSize*blockSize, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&BlockB[0][0], blockSize* blockSize, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		}
		for (int i = 0; i < sqrtAB; i++) {
			addMult(BlockA, BlockB, blockSize, BlockC);
			if (pID == closeToRandomVariable) {
				/*cout << "I multiplied a*b=c " << endl;
				printMatrix(BlockA, blockSize, blockSize, "BlockA");
				printMatrix(BlockB, blockSize, blockSize, "BlockB");
				printMatrix(BlockC, blockSize, blockSize, "BlockC");*/
			}
			int row_Dest, col_Dest;
			row_Dest = (row_rank - 1);
			row_Dest = (row_Dest + sqrtAB) % sqrtAB;
			col_Dest = (col_rank - 1);
			col_Dest = (col_Dest + sqrtAB) % sqrtAB;
			//std::copy(&BlockA[0][0],&BlockA[0][0]+blockSize*blockSize,&a2[0][0]);
			//std::copy(&BlockB[0][0],&BlockB[0][0]+blockSize*blockSize,&b2[0][0]);
			MPI_Sendrecv(&BlockA[0][0], blockSize*blockSize, MPI_DOUBLE, row_Dest, 0, &a2[0][0], blockSize*blockSize, MPI_DOUBLE, (row_rank + 1) % sqrtAB, 0, rowCom, &statusForA);
			MPI_Sendrecv(&BlockB[0][0], blockSize*blockSize, MPI_DOUBLE, col_Dest, 0, &b2[0][0], blockSize*blockSize, MPI_DOUBLE, (col_rank + 1) % sqrtAB, 0, colCom, &statusForB);
			std::copy(&a2[0][0], &a2[0][0] + blockSize*blockSize, &BlockA[0][0]);
			std::copy(&b2[0][0], &b2[0][0] + blockSize*blockSize, &BlockB[0][0]);
		}
		if (pID != 0) {
			MPI_Send(&BlockC[0][0], blockSize*blockSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
//----------------------------Recv Blocks------------------------//
	if (pID == 0) {
		double** C = alloc_2d_int(totalMSize, totalMSize);
		for (int i = 0; i < sqrtAB; i++) {
			for (int j = 0; j < sqrtAB; j++) {
				int destination = j + (i * sqrtAB);
				if (destination != 0) {
					MPI_Recv(&BlockC[0][0], blockSize*blockSize, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD, &status);
				}
				blockIntoMat(BlockC, C, i*blockSize, j*blockSize, blockSize, blockSize);
			}
		}
		std::chrono::system_clock::time_point endTime = std::chrono::system_clock::now();
		std::chrono::microseconds microRunTime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
		double runTime = microRunTime.count() / 1000000.0;
		std::cout << std::endl;
		for (int r = 0; r < totalMSize; r++) {
			for (int c = 0; c < totalMSize; c++) {
				cout << C[r][c] << " ";
			}
			cout << endl;
		}	
		std::cout << "Wall clock time = " << runTime << " seconds." << std::endl << std::flush;
	}
	MPI_Finalize();	
    return 0;
}
