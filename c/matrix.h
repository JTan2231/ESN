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
#include <complex.h>
#include "generation.h"
#include "map.h"

typedef struct {
    double** array;
    int rows;
    int cols;
} Matrix;

typedef struct {
    double* array;
    int size;
} Vector;

// TODO: Organizing this file
// TODO: Tikhonov Regularization

//----------------------------------------\\
// Initialization                         \\
//----------------------------------------\\

// initialize matrix of zeros of size (r, c)
void initMat(Matrix* mat, int r, int c) {
    mat->array = calloc(r, sizeof(*(mat->array)));
    for (int i = 0; i < r; i++)
        mat->array[i] = calloc(c, sizeof(*(mat->array[i])));
    
    mat->rows = r;
    mat->cols = c;
}

// initialize empty vector of given size
void initVec(Vector* vec, int size) {
    vec->array = calloc(size, sizeof(vec->array));
    vec->size = size;
}

// initialize identity matrix of size (r, c)
void initIdent(Matrix* mat, int r, int c) {
    mat->array = calloc(r, sizeof(*(mat->array)));
    for (int i = 0; i < r; i++)
        mat->array[i] = calloc(c, sizeof(*(mat->array[i])));
    
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++)
        mat->array[i][i] = 1;
}

// deallocate the given matrix
void cleanMat(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++)
        free(mat->array[i]);
        
    free(mat->array);
}

// deallocate the given vector
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

// initialize vector with random floating point values
void initVecRandom(Vector* vec, int size) {
    vec->array = calloc(size, sizeof(vec->array));
    vec->size = size;
    
    for (int i = 0; i < size; i++) 
        vec->array[i] = randomDouble();
}

// initialize matrix with random integers between 0 - range
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
// uses Sparse as a list of columns
// each nested Column acting as a unique column
// with a list of rows in which there is a value
// for each row where there is a value in the sparse matrix
// eliminates unnecessary storage of zeros
void initSparse(Sparse* mat, int r, int c, double density) {
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
        int in = sparseIn(mat, column);
        if (in+1 != mat->arraySize) {
            // if the generated column number exists, check for row
            int row = randRange(r);
            int in2 = columnInFirst(&(mat->array[in]), row);
            if (in2 == 0) {
                // row doesn't exist, add the value
                columnAdd(&(mat->array[in]), row, marsagliaPolar());
                count++;
            }
        }
        else {
            // make and add a new column
            Column p;
            initColumn(&p, column);
            columnAdd(&p, randRange(r), marsagliaPolar());
            sparseAdd(mat, &p);
            count++;
        }
    }
}

// initialize matrix using normally distributed random values
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

// initialize vector using normally distributed random values
void initVecRandomNormal(Vector* vec, int size) {
    vec->array = calloc(size, sizeof(vec->array));
    vec->size = size;
    
    printf("check\n");
    for (int i = 0; i < size; i++)
        vec->array[i] = marsagliaPolar();
}

// initialize matT to be the transpose of mat
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

void initClone(Matrix* mat, Matrix* clone) {
    clone->array = calloc(mat->rows, sizeof(*(mat->array)));
    for (int i = 0; i < mat->rows; i++)
        clone->array[i] = calloc(mat->cols, sizeof(*(mat->array[i])));
    
    clone->rows = mat->rows;
    clone->cols = mat->cols;

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            clone->array[i][j] = mat->array[i][j];
    }
}

void initVecRow(Matrix* mat, int row, Vector* vec) {
    vec->array = calloc(mat->rows, sizeof(vec->array));
    vec->size = mat->rows;
    
    for (int i = 0; i < mat->rows; i++)
        vec->array[i] = mat->array[row][i];
}

// initializes a vector from a matrix column
void initVecCol(Matrix* mat, int col, Vector* vec) {
    vec->array = calloc(mat->rows, sizeof(vec->array));
    vec->size = mat->rows;
    
    for (int i = 0; i < mat->rows; i++)
        vec->array[i] = mat->array[i][col];
}

void initSparseToMat(Sparse* sparse, Matrix* mat) {
    mat->array = calloc(sparse->rows, sizeof(mat->array));
    for (int i = 0; i < sparse->rows; i++)
        mat->array[i] = calloc(sparse->cols, sizeof(*(mat->array[i])));

    mat->rows = sparse->rows;
    mat->cols = sparse->cols;

    for (int j = 0; j < sparse->arraySize; j++) {
        for (int i = 0; i < sparse->array[j].arraySize; i++) {
            int matRow = sparse->array[j].array[i].first;
            int matCol = sparse->array[j].value;
            mat->array[matRow][matCol] = sparse->array[j].array[i].second;
        }
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

void printVec(Vector* vec) {
    printf("[ \n");
    for (int i = 0; i < vec->size; i++)
        printf("\t%.16lf \n", vec->array[i]);
    printf("]\n");
}

//--------------------------------------------------------\\
// Arithmetic Operations                                  \\
//--------------------------------------------------------\\

void vecSub(Vector* v1, Vector* v2) {
    assert(v1->size == v2->size);

    for (int i = 0; i < v1->size; i++)
        v1->array[i] -= v2->array[i];
}

void vecAdd(Vector* v1, Vector* v2) {
    assert(v1->size == v2->size);
    
    for (int i = 0; i < v1->size; i++)
        v1->array[i] += v2->array[i];
}

// multiply all elements of given vector
// by given scalar
void scalarVec(Vector* vec, double scalar) {
    for (int i = 0; i < vec->size; i++)
        vec->array[i] *= scalar;
}

// vec = vec - scalar
void scalarVecSub(Vector* vec, double scalar) {
    for (int i = 0; i < vec->size; i++)
        vec->array[i] -= scalar;
}

// vec = vec / scalar
void scalarVecDivOut(Vector* vec, double scalar, Vector* out) {
    for (int i = 0; i < vec->size; i++)
        out->array[i] = vec->array[i] / scalar;
}

// component function
// for use in matDot
// dot product of m1 row and m2 column
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

double matVecDotPartial(Matrix* mat, int row, Vector* vec) {
    double output = 0;
    for (int i = 0; i < vec->size; i++)
        output += mat->array[row][i] * vec->array[i];
    
    return output;
}

double matVecDotPartialCol(Matrix* mat, int col, Vector* vec) {
    double output = 0;
    for (int i = 0; i < vec->size; i++)
        output += mat->array[i][col] * vec->array[i];

    return output;
}

void matVecDot(Matrix* mat, Vector* vec, Vector* out) {
    assert(mat->cols == vec->size);
    assert(mat->rows == out->size);

    for (int i = 0; i < out->size; i++)
        out->array[i] = matVecDotPartial(mat, i, vec);
}

// component function
// for use in sparseDotSecond
// computes the dot product of 
// a dense matrix row 
// and sparse matrix column
double sparseDotPartial(Matrix* mat, int row, Sparse* sparse, int col) {
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
void sparseVecDot(Sparse* sparse, Vector* vec, Vector* out) {
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
void sparseDotSecond(Matrix* mat, Sparse* sparse, Matrix* out) {
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
double vecDot(Vector* col, Vector* row) {
    double out = 0;
    for (int i = 0; i < col->size; i++)
        out += col->array[i] * row->array[i];

    return out;
}

//--------------------------------------------------------\\
// Manipulation                                           \\
//--------------------------------------------------------\\

// sets the given spot to a certain value
// unused for some reason
void set(Matrix* mat, double value, int row, int col) {
    mat->array[row][col] = value;
}

// assign matrix column from vector
void assignMatColVec(Matrix* mat, int col, Vector* vec) {
    assert(mat->rows == vec->size);
    
    for (int i = 0; i < mat->rows; i++)
        mat->array[i][col] = vec->array[i];
}

//--------------------------------------------------------\\
// Linear Algebra                                         \\
//--------------------------------------------------------\\

// normalizes the given vector
/*void normalize(Vector* vec, double magnitude) {
    for (int i = 0; i < vec->size; i++)
        vec->array[i] /= magnitude;
}

// normalizes the given vector
// and outputs the normalized values
// in the given output vector
void normalizeOut(Vector* vec, Vector* output, double magnitude) {
    assert(vec->size == output->size);

    assert(magnitude > 0);
    
    for (int i = 0; i < vec->size; i++) {
        double assigned = vec->array[i] / magnitude;
        assert(!isnan(assigned));
        output->array[i] = assigned;
    }
}

// calculates the magnitude of the given vector
double magnitude(Vector* vec) {
    double output = 0;
    for (int i = 0; i < vec->size; i++)
        output += vec->array[i] * vec->array[i];

    return sqrt(output);
}

// uses power iteration to find the eigenVector
// for a sparse matrix
// see https://en.wikipedia.org/wiki/Power_iteration
void eigenVectorSparse(Sparse* mat, Vector* eigenVec, int iterations) {
    Vector dot;
    for (int i = 0; i < iterations; i++) {
        initVec(&dot, eigenVec->size);
        sparseVecDot(mat, eigenVec, &dot);
        normalizeOut(&dot, eigenVec, magnitude(&dot));
        cleanVec(&dot);
    }
}

// uses power iteration to find the eigenVector
// for a dense matrix
// see https://en.wikipedia.org/wiki/Power_iteration
void eigenVectorDense(Matrix* mat, Vector* eigenVec, int iterations) {
    Vector dot;
    for (int i = 0; i < iterations; i++) {
        initVec(&dot, eigenVec->size);
        matVecDot(mat, eigenVec, &dot);
        normalizeOut(&dot, eigenVec, magnitude(&dot));
        cleanVec(&dot);
    }
}

// uses the Rayleigh Quotient to get the spectral radius
// of the given eigenvector
// presumably from power iteration
// see https://en.wikipedia.org/wiki/Power_iteration
double rayleighQuotient(Sparse* sparse, Vector* vec) {
    Vector sparseDot;
    initVec(&sparseDot, vec->size);
    sparseVecDot(sparse, vec, &sparseDot);

    double num = vecDot(vec, &sparseDot);
    double den = vecDot(vec, vec);

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
}*/

#endif
