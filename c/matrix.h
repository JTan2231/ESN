/*
 * Author: Joey Tan
 * Date Created: 2-18-20
 * Last Edit: 3-4-20, Joey Tan
 */

#ifndef MATRIX
#define MATRIX
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "generation.h"
#include "map.h"
#include "complex.h"

typedef struct {
    double** array;
    int rows;
    int cols;
} Matrix;

typedef struct {
    double* array;
    int size;
} Vector;

typedef struct {
    Complex* ptr;
    int rows;
    int cols;
} ComplexMatrix;

// TODO: Organizing this file
// TODO: Tikhonov Regularization

//----------------------------------------\\
// Initialization                         \\
//----------------------------------------\\

// initialize matrix of zeros
void initMat(Matrix* mat, int r, int c) {
    mat->array = calloc(r, sizeof(*(mat->array)));
    for (int i = 0; i < r; i++)
        mat->array[i] = calloc(c, sizeof(*(mat->array[i])));
    
    mat->rows = r;
    mat->cols = c;
}

void initVec(Vector* vec, int size) {
    vec->array = calloc(size, sizeof(vec->array));
    vec->size = size;
}

void initIdent(Matrix* mat, int r, int c) {
    mat->array = calloc(r, sizeof(*(mat->array)));
    for (int i = 0; i < r; i++)
        mat->array[i] = calloc(c, sizeof(*(mat->array[i])));
    
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++)
        mat->array[i][i] = 1;
}

void cleanMat(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++)
        free(mat->array[i]);
        
    free(mat->array);
}

void cleanVec(Vector* vec) {
    free(vec->array);
}

// initialize matrix with random floating point values
void initRandom(Matrix* mat, int r, int c) {
    mat->array = calloc(r, sizeof(*(mat->array)));
    for (int i = 0; i < r; i++)
        mat->array[i] = calloc(c, sizeof(*(mat->array[i])));
    
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->array[i][j] = randomDouble();
    }
}

void initVecRandom(Vector* vec, int size) {
    vec->array = calloc(size, sizeof(vec->array));
    vec->size = size;
    
    for (int i = 0; i < size; i++) 
        vec->array[i] = randomDouble();
}

// initialize a matrix with random integers between 0 - range
void initRandomRange(Matrix* mat, int r, int c, int range) {
    mat->array = calloc(r, sizeof(*(mat->array)));
    for (int i = 0; i < r; i++)
        mat->array[i] = calloc(c, sizeof(*(mat->array[i])));
    
    mat->rows = r;
    mat->cols = c;
    
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->array[i][j] = randRange(range);
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
    mat->array = calloc(r, sizeof(*(mat->array)));
    for (int i = 0; i < r; i++)
        mat->array[i] = calloc(c, sizeof(*(mat->array[i])));
    
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->array[i][j] = marsagliaPolar();
    }
}

void initVecRandomNormal(Vector* vec, int size) {
    vec->array = calloc(size, sizeof(vec->array));
    vec->size = size;
    
    for (int i = 0; i < size; i++)
        vec->array[i] = marsagliaPolar();
}

// initializes matT to be the transpose of mat
void initTranspose(Matrix* mat, Matrix* matT) {
    matT->array = calloc(mat->rows, sizeof(*(mat->array)));
    for (int i = 0; i < mat->rows; i++)
        matT->array[i] = calloc(mat->cols, sizeof(*(mat->array[i])));
    
    matT->rows = mat->rows;
    matT->cols = mat->cols;

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            matT->array[j][i] = mat->array[i][j];
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
            printf("\t%.16lf \n", mat->array[i][j]);
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
        output += m1->array[row][i] * m2->array[i][col];

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
            m3->array[i][j] = matDotPartial(m1, i, m2, j);
    }
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
        double matrixValue = mat->array[row][sparseColIndex];
        double sparseValue = sparse->array[col].array[i].second;
        
        output += matrixValue * sparseValue;
    }

    return output;
}

// matrix multiplication == sparseMatrix * columnVector
// This one is specifically for the above expression
// it is not generalized for all matrices
void sparseDotFirst(ParentArray* sparse, Vector* vec, Vector* out) {
    assert(out->size == vec->size || out->size == sparse->cols || vec->size == sparse->rows);

    for (int i = 0; i < sparse->arraySize; i++) {
        for (int j = 0; j < sparse->array[i].arraySize; j++) {
            int sparseCol = sparse->array[i].value;
            double sparseValue = sparse->array[i].array[j].second;
            double matrixValue = vec->array[sparse->array[i].value];

            out->array[sparseCol] +=  sparseValue * matrixValue;
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
            out->array[i][sparse->array[j].value] = sparseDotPartial(mat, i, sparse, j);
    }
}

// vector dot == columnVector dot rowVector
// dot product 
// of matrix representations 
// of vectors
double colRowDot(Vector* col, Vector* row) {
    double out = 0;
    for (int i = 0; i < col->size; i++)
        out += col->array[i] * row->array[i];

    return out;
}

//--------------------------------------------------------\\
// Manipulation                                           \\
//--------------------------------------------------------\\

void set(Matrix* mat, double value, int row, int col) {
    mat->array[row][col] = value;
}

//--------------------------------------------------------\\
// Linear Algebra                                         \\
//--------------------------------------------------------\\

// only meant to be used on vectors (1-D matrices)
// normalizes the given vector
void normalize(Vector* vec, double magnitude) {
    assert(vec->size == 1);

    for (int i = 0; i < vec->size; i++)
        vec->array[i] /= magnitude;
}

// only meant to be used on vectors (1-D matrices)
// normalizes the given vector
// and outputs the normalized values
// in the given output vector
void normalizeOut(Vector* vec, Vector* output, double magnitude) {
    assert(vec->size == 1);
    assert(output->size == 1);
    assert(vec->size == output->size);

    assert(magnitude > 0);
    
    for (int i = 0; i < vec->size; i++) {
        double assigned = vec->array[i] / magnitude;
        assert(!isnan(assigned));
        output->array[i] = assigned;
    }
}

// only meant to be used on vectors (1-D matrices)
// calculates the magnitude of the given vector
double magnitude(Vector* vec) {
    assert(vec->size == 1);

    double output = 0;
    for (int i = 0; i < vec->size; i++)
        output += vec->array[i] * vec->array[i];

    return sqrt(output);
}

// uses power iteration to find the eigenVector
// for the sparse matrix
// see https://en.wikipedia.org/wiki/Power_iteration
void eigenVector(ParentArray* mat, Vector* eigenVec, int iterations) {
    Vector dot;
    for (int i = 0; i < iterations; i++) {
        initVec(&dot, eigenVec->size);
        sparseDotFirst(mat, eigenVec, &dot);
        normalizeOut(&dot, eigenVec, magnitude(&dot));
        cleanVec(&dot);
    }
}

// uses the Rayleigh Quotient to get the spectral radius
// of the given eigenvector
// presumably from power iteration
// see https://en.wikipedia.org/wiki/Power_iteration
double rayleighQuotient(ParentArray* sparse, Vector* vec) {
    Vector sparseDot;
    initVec(&sparseDot, vec->size);
    sparseDotFirst(sparse, vec, &sparseDot);

    double num = colRowDot(vec, &sparseDot);
    double den = colRowDot(vec, vec);

    return num / den;
}

// computes LU decomposition
// of given matrix (Crout style)
// see https://en.wikipedia.org/wiki/Crout_matrix_decomposition
void crout(Matrix* mat, Matrix* lower, Matrix* upper) {
    assert(mat->rows == mat->cols);

    for (int i = 0; i < mat->rows; i++)
        upper->array[i][i] = 1;

    // for loop hell
    for (int i = 0; i < mat->rows; i++) {
        for (int j = i; j < mat->rows; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++)
                sum += lower->array[j][k] * upper->array[k][i];

            lower->array[j][i] = mat->array[j][i] - sum;
        }

        for (int j = i; j < mat->rows; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++)
                sum += lower->array[i][k] * upper->array[k][j];

            upper->array[i][j] = (mat->array[i][j] - sum) / lower->array[i][i];
        }
    }
}

// computes the inverse of the given matrix
// using front- and back-substitution
// see: 
// - https://algowiki-project.org/en/Backward_substitution
// - https://algowiki-project.org/en/Forward_substitution
void inverse(Matrix* mat, Matrix* lower, Matrix* upper, Matrix* inverse) {
    Matrix eye;
    initIdent(&eye, mat->rows, mat->cols);

    Matrix d;
    initMat(&d, mat->rows, mat->cols);

    for (int column = 0; column < d.cols; column++) {
        for (int i = 0; i < d.rows; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += lower->array[i][j] * d.array[j][column];

            d.array[i][column] = (eye.array[i][column]-sum) / lower->array[i][i];
        }
    }

    // solve for x (the inverse of mat)
    for (int column = 0; column < inverse->cols; column++) {
        for (int i = inverse->rows-1; i >= 0; i--) {
            double sum = 0;
            for (int j = i; j < inverse->rows; j++)
                sum += upper->array[i][j] * inverse->array[j][column];

            inverse->array[i][column] = d.array[i][column] - sum;
        }
    }
}


#endif
