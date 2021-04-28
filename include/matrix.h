/*
 * Author: Joey Tan
 * Date Created: 2-18-20
 * Last Edit: 3-13-20, Joey Tan
 */

#ifndef MATRIX
#define MATRIX
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <complex.h>
#include "generation.h"
#include "sparse.h"

typedef struct {
    double** array;
    int rows;
    int cols;
    int size;
} Matrix;

typedef struct {
    double* array;
    int size;
} Vector;

// TODO: Organizing this file
// TODO: Documentation on what these functions do

//---------------------------------------------\\
// Matrix Initialization                       \\
//                                             \\
// Functions to initialize a matrix.           \\
// General format:                             \\
// -- init___(Matrix* mat, int r, int c);      \\
// where:                                      \\
// -- mat: pointer to matrix being initialized \\
// -- -- NOTE: initializing a matrix using an  \\
//             already calloc'ed pointer       \\
//             leaves the previous matrix      \\
//             data unaccessable               \\
//                                             \\
// -- r: desired number of rows                \\
// -- c: desired number of columns             \\
//---------------------------------------------\\

// initializes an empty matrix of zeros
void initMat(Matrix* mat, int r, int c);
// initializes a matrix of ones
void initOnes(Matrix* mat, int r, int c);
// initializes an identity matrix
void initIdent(Matrix* mat, int r, int c);

// initializes Sparse struct from sparse.h
void initSparse(Sparse* mat, int r, int c, double density);

// evenly distributed real numbers instead of zeros
// -1 <= x <= 1
void initRandom(Matrix* mat, int r, int c);

// normally distributed real numbers
// mu == 0, variance == 1
void initRandomNormal(Matrix* mat, int r, int c);

// initializes a matrix (matT) to be the transpose
// of an already-initialized matrix (mat)
void initTranspose(Matrix* mat, Matrix* matT);

// initializes a matrix (clone) to be a clone 
// of an already-initialized matrix (mat)
void initClone(Matrix* mat, Matrix* clone);

// initializes a matrix (mat) from an already-initialized
// sparse matrix (sparse)
void initSparseToMat(Sparse* sparse, Matrix* mat);

// These functions clean and initialize their respective structs
void reinitMat(Matrix* mat, int r, int c);
void reinitSparse(Sparse* sparse, int r, int c, double d);

//----------------------------------------\\
// Vector Initialization                  \\
//                                        \\
// Basically the same as above            \\
//----------------------------------------\\

void initVec(Vector* vec, int size);
void initVecRandom(Vector* vec, int size);
void initVecRandomNormal(Vector* vec, int size);

void initVecClone(Vector* vec, Vector* clone);

// initialize a vector using values from the row
// of a matrix
void initVecRow(Matrix* mat, int row, Vector* vec);
// above but column
void initVecCol(Matrix* mat, int col, Vector* vec);
// above but only part of the column
void initVecColPartial(Matrix* mat, int rowBegin, int rowEnd, int col, Vector* vec);

//----------------------------------------\\
// Matrix and Vector Operations           \\
//                                        \\
// Matrix and vector arithmetic.          \\
//----------------------------------------\\

// Vector

// v1 -= v2
void vecSub(Vector* v1, Vector* v2);
// out = v1 -v2
void vecSubOut(Vector* v1, Vector* v2, Vector* out);

// v1 += v2
void vecAdd(Vector* v1, Vector* v2);
// out = v1 + v2
void vecAddOut(Vector* v1, Vector* v2, Vector* out);

// vec *= scalar
void scalarVec(Vector* vec, double scalar);
// out = vec * scalar
void scalarVecOut(Vector* vec, double scalar, Vector* out);

// vec -= scalar
void scalarVecSub(Vector* vec, double scalar);

// vec /= scalar
void scalarVecDiv(Vector* vec, double scalar);
// out = vec / scalar
void scalarVecDivOut(Vector* vec, double scalar, Vector* out);

// vector dot product
double vecDot(Vector* col, Vector* row);

// vector dot products with a matrix row or column
double matVecDotPartialRow(Matrix* mat, int row, Vector* vec);
double matVecDotPartialCol(Matrix* mat, int col, Vector* vec);

// out = mat * vec
void matVecDot(Matrix* mat, Vector* vec, Vector* out);
// out = vec * mat
void vecMatDot(Vector* vec, Matrix* mat, Vector* out);

// the above but with a sparse matrix
void sparseVecDot(Sparse* sparse, Vector* vec, Vector* out);
void vecSparseDot(Vector* vec, Sparse* sparse, Vector* out);

// Matrix

// out = m1 + m2
void matAdd(Matrix* m1, Matrix* m2, Matrix* out);
// m3 = m1 * m2
void matDot(Matrix* m1, Matrix* m2, Matrix* m3);

// mat ?= scalar, where ? is the operation (+, -, etc.)
void scalarMatDiv(Matrix* mat, double scalar);
void scalarMatSub(Matrix* mat, double scalar);
void scalarMatMult(Matrix* mat, double scalar);
void scalarSparseMult(Sparse* mat, double scalar);
void scalarSparseDiv(Sparse* mat, double scalar);

// matDot but for sparse matrices
void sparseDotFirst(Sparse* sparse, Matrix* mat, Matrix* out);
void sparseDotSecond(Matrix* mat, Sparse* sparse, Matrix* out);

// vector dot product with matrix rows/columns
double matDotPartial(Matrix* m1, int row, Matrix* m2, int col);
double sparseDotPartial(Matrix* mat, int row, Sparse* sparse, int col);

//----------------------------------------\\
// Manipulation/Miscellaneous             \\
//----------------------------------------\\

// Cleanup
void cleanMat(Matrix* mat);
void cleanVec(Vector* vec);

// I/O
void printMat(Matrix* mat);
void printVec(Vector* vec);

// Checks
int checkIdent(Matrix* mat);
int checkConvert(Matrix* mat, Vector* vec);

// Conversion

// TODO: rename these for clarity
// converts a Matrix struct into a Vector struct
// NOTE: original matrix/vector remains after conversion
void convertMatToVec(Matrix* mat, Vector* vec);
// and vice versa
void convertVecToMat(Vector* vec, Matrix* mat);

// Assignment

// overwrites a vector with all the elements of a sparse matrix row
void assignVecSparseRow(Vector* vec, int row, Sparse* sparse);
// overwrites matrix column with all vector elements
void assignMatColVec(Matrix* mat, int col, Vector* vec);
// the first but with a Matrix struct and a column
void assignVecMatCol(Vector* vec, int col, Matrix* mat);

// assigns all elements of a vector/matrix to 0
void zeroMat(Matrix* mat);
void zeroVec(Vector* vec);

// assigns all elements of m1/v2 to all elements of m2/v2
void cloneMat(Matrix* m1, Matrix* m2);
void cloneVec(Vector* v1, Vector* v2);

// appends a matrix onto another matrix, row-wise
void appendMatRow(Matrix* m1, Matrix* m2);
// appends a column vector to a matrix
void appendMatColVec(Matrix* mat, Vector* vec);

// reduces size of matrix
// values outside new dimension are discarded
void shrinkMat(Matrix* mat, int rows, int cols);

// increases size of matrix
// new values are 0
void growMat(Matrix* mat, int rows, int cols);

// creates a Matrix struct object from a Sparse struct object
void sparseToMat(Sparse* sparse, Matrix* mat);

#endif
