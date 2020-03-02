/*
 * Author: Joey Tan
 * Date Created: 2-18-20
 * Last Edit: 2-23-20, Joey Tan
 */

#ifndef MATRIX
#define MATRIX
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "generation.h"
#include "map.h"

typedef struct {
    double* ptr;
    int rows;
    int cols;
} Matrix;

// TODO: Organizing this file
// TODO: Matrix inverse
// TODO: Tikhonov Regularization
// TODO: Vector structure and functions
// TODO: Matrix pointer to pointer pointer for easier indexing

//----------------------------------------\\
// Initialization                         \\
//----------------------------------------\\

// initialize matrix of zeros
void initMat(Matrix* mat, int r, int c) {
    mat->ptr = calloc(r * c, sizeof mat->ptr);
    mat->rows = r;
    mat->cols = c;
}

void cleanMat(Matrix* mat) {
    free(mat->ptr);
}

// initialize matrix with random floating point values
void initRandom(Matrix* mat, int r, int c) {
    mat->ptr = calloc(r * c, sizeof mat->ptr);
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->ptr[i*mat->cols+j] = randomDouble();
    }
}

// initialize a matrix with random integers between 0 - range
void initRandomRange(Matrix* mat, int r, int c, int range) {
    mat->ptr = calloc(r * c, sizeof mat->ptr);
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->ptr[i*mat->cols+j] = randRange(range);
    }
}

// uses normally distributed values
// uses ParentArray as a list of columns
// each nested Parent acting as a unique column
// with a list of rows in which there is a value
// for each row where there is a value in the sparse matrix
// eliminates unnecessary storage of zeros
void initSparse(ParentArray* mat, int r, int c, double density) {
    mat->array = malloc(sizeof mat->array);
    mat->arraySize = 1;
    mat->density = density;
    mat->rows = r;
    mat->cols = c;
    
    int count = 0;
    int total = (int)(r*c*density);

    // density is just a percentage
    // of the max available
    // to be filled
    while (count < total) {
        int column = randRange(c);
        int in = parentArrayIn(mat, column);
        if (in+1 != mat->arraySize) {
            // if the generated column number exists, check for row
            int row = randRange(r);
            int in2 = parentInFirst(&(mat->array[in]), row);
            if (in2 == 0) {
                parentAdd(&(mat->array[in]), row, marsagliaPolar());
                count++;
            }
        }
        else {
            // make and add a new column
            Parent p;
            initParent(&p, column);
            parentAdd(&p, randRange(r), marsagliaPolar());
            parentArrayAdd(mat, &p);
            count++;
        }
    }
}

// initializes a matrix using normally distributed random values
void initRandomNormal(Matrix* mat, int r, int c) {
    mat->ptr = calloc(r * c, sizeof mat->ptr);
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->ptr[i*mat->cols+j] = marsagliaPolar();
    }
}

// initializes matT to be the transpose of mat
void initTranspose(Matrix* mat, Matrix* matT) {
    matT->ptr = calloc(mat->rows * mat->cols, sizeof matT->ptr);
    matT->rows = mat->cols;
    matT->cols = mat->rows;

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            matT->ptr[j*matT->cols+i] = mat->ptr[i*mat->cols+j];
    }
}

//--------------------------------------------------------\\
// I/O                                                    \\
//--------------------------------------------------------\\

// prints the given matrix
void printMat(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        printf("[ \n");
        for (int j = 0; j < mat->cols; j++)
            printf("\t%.16lf \n", mat->ptr[i*mat->cols+j]);
        printf("]\n");
    }
}



//--------------------------------------------------------\\
// Arithmetic Operations                                  \\
//--------------------------------------------------------\\

// component function
// for use in matDot
double matDotPartial(Matrix* m1, int row, Matrix* m2, int col) {
    double output = 0;
    for (int i = 0; i < m2->rows; i++)
        output += m1->ptr[row*m1->cols+i] * m2->ptr[i*m2->cols+col];

    return output;
}

// general matrix multiplication function
// not for use with vectors
void matDot(Matrix* m1, Matrix* m2, Matrix* m3) {
    assert(m1->cols == m2->rows);
    assert(m3->rows == m1->rows);
    assert(m3->cols == m2->cols);

    for (int i = 0; i < m1->rows; i++) {
        for (int j = 0; j < m2->cols; j++)
            m3->ptr[i*m3->cols+j] = matDotPartial(m1, i, m2, j);
    }
}

// matrix by vector dot product
void matVecDot(Matrix* mat, Matrix* vec, Matrix* outVec) {
    assert(mat->rows == vec->cols);
    assert(vec->cols == outVec->cols);

    for (int i = 0; i < mat->rows; i++)
        outVec->ptr[i] = mat->ptr[i*mat->cols] * vec->ptr[i];
}

// component function
// for use in sparseDotSecond
// computes the dot product of 
// a dense matrix row 
// and sparse matrix column
double sparseDotPartial(Matrix* mat, int row, ParentArray* sparse, int col) {
    double output = 0;
    for (int i = 0; i < sparse->array[col].arraySize; i++) {
        int sparseColIndex = sparse->array[col].array[i].first;
        double matrixValue = mat->ptr[row*mat->cols+sparseColIndex];
        double sparseValue = sparse->array[col].array[i].second;
        
        output += matrixValue * sparseValue;
    }

    return output;
}

// matrix multiplication == sparseMatrix * columnVector
// This one is specifically for the above expression
// it is not generalized for all matrices
void sparseDotFirst(ParentArray* sparse, Matrix* mat, Matrix* out) {
    assert(out->rows == mat->rows || out->cols == sparse->cols || mat->cols == sparse->rows);

    for (int i = 0; i < sparse->arraySize; i++) {
        for (int j = 0; j < sparse->array[i].arraySize; j++) {
            int sparseCol = sparse->array[i].value;
            double sparseValue = sparse->array[i].array[j].second;
            double matrixValue = mat->ptr[sparse->array[i].value];

            out->ptr[sparseCol] +=  sparseValue * matrixValue;
        }
    }
}

// matrix multiplication == denseMatrix * sparseMatrix
void sparseDotSecond(Matrix* mat, ParentArray* sparse, Matrix* out) {
    assert(mat->cols == sparse->rows);
    assert(out->rows == mat->rows);
    assert(out->cols == sparse->cols);

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < sparse->arraySize; j++)
            out->ptr[i*out->cols+sparse->array[j].value] = sparseDotPartial(mat, i, sparse, j);
    }
}

// vector dot == columnVector dot rowVector
// dot product 
// of matrix representations 
// of vectors
double colRowDot(Matrix* col, Matrix* row) {
    double out = 0;
    for (int i = 0; i < col->rows; i++)
        out += col->ptr[i] * row->ptr[i];

    return out;
}

//--------------------------------------------------------\\
// Manipulation                                           \\
//--------------------------------------------------------\\

void set(Matrix* mat, double value, int row, int col) {
    mat->ptr[row*mat->cols+col] = value;
}

//--------------------------------------------------------\\
// Linear Algebra                                         \\
//--------------------------------------------------------\\

// only meant to be used on vectors (1-D matrices)
// normalizes the given vector
void normalize(Matrix* vec, double magnitude) {
    assert(vec->rows == 1);

    for (int i = 0; i < vec->cols; i++)
        vec->ptr[i] /= magnitude;
}

// only meant to be used on vectors (1-D matrices)
// normalizes the given vector
// and outputs the normalized values
// in the given output vector
void normalizeOut(Matrix* vec, Matrix* output, double magnitude) {
    assert(vec->rows == 1);
    assert(output->rows == 1);
    assert(vec->cols == output->cols);

    assert(magnitude > 0);
    
    for (int i = 0; i < vec->cols; i++) {
        double assigned = vec->ptr[i] / magnitude;
        assert(!isnan(assigned));
        output->ptr[i] = assigned;
    }
}

// only meant to be used on vectors (1-D matrices)
// calculates the magnitude of the given vector
double magnitude(Matrix* vec) {
    assert(vec->rows == 1);

    double output = 0;
    for (int i = 0; i < vec->cols; i++)
        output += vec->ptr[i] * vec->ptr[i];

    return sqrtf(output);
}

// uses power iteration to find the eigenVector
// for the sparse matrix
// see https://en.wikipedia.org/wiki/Power_iteration
void eigenVector(ParentArray* mat, Matrix* eigenVec, int iterations) {
    Matrix dot;
    for (int i = 0; i < iterations; i++) {
        initMat(&dot, 1, eigenVec->cols);
        sparseDotFirst(mat, eigenVec, &dot);
        normalizeOut(&dot, eigenVec, magnitude(&dot));
        cleanMat(&dot);
    }
}

// uses the Rayleigh Quotient to get the spectral radius
// of the given eigenvector
// presumably from power iteration
// see https://en.wikipedia.org/wiki/Power_iteration
double rayleighQuotient(ParentArray* sparse, Matrix* vec) {
    Matrix vecT;
    initTranspose(vec, &vecT);
    
    Matrix sparseDot;
    initMat(&sparseDot, vec->rows, vec->cols);
    sparseDotFirst(sparse, vec, &sparseDot);

    double num = colRowDot(&vecT, &sparseDot);
    double den = colRowDot(&vecT, vec);

    return num / den;
}

#endif
