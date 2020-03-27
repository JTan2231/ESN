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
// TODO: EIGENVALUES
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

// initialize vector with random floating point values
void initVecRandom(Vector* vec, int size) {
    vec->array = calloc(size, sizeof(vec->array));
    vec->size = size;
    
    for (int i = 0; i < size; i++) 
        vec->array[i] = randomDouble();
}

// deallocate the given matrix
void cleanMat(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++)
        free(mat->array[i]);
        
    free(mat->array);
}

// set all elements in given matrix to zero
void zeroMat(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            mat->array[i][j] = 0;
    }
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
                columnAdd(&(mat->array[in]), row, randomDouble());
                count++;
            }
        }
        else {
            // make and add a new column
            Column p;
            initColumn(&p, column);
            columnAdd(&p, randRange(r), randomDouble());
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
    
    for (int i = 0; i < size; i++)
        vec->array[i] = marsagliaPolar();
}

// initialize matT to be the transpose of mat
void initTranspose(Matrix* mat, Matrix* matT) {
    matT->array = calloc(mat->rows, sizeof(*(mat->array)));
    for (int i = 0; i < mat->rows; i++)
        matT->array[i] = calloc(mat->cols, sizeof(*(mat->array[i])));
    
    matT->rows = mat->cols;
    matT->cols = mat->rows;

    for (int i = 0; i < matT->rows; i++) {
        for (int j = 0; j < matT->cols; j++)
            matT->array[i][j] = mat->array[j][i];
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

void reinitMat(Matrix* mat, int r, int c) {
    cleanMat(mat);
    initMat(mat, r, c);
}

void reinitSparse(Sparse* sparse, int r, int c, double d) {
    //printf("Reinitialization check. density: %lf\n", d);
    cleanSparse(sparse);
    initSparse(sparse, r, c, d);
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

void convertMatToVec(Matrix* mat, Vector* vec) {
    int a = mat->rows == vec->size;
    int b = mat->cols == 1;
    int c = mat->rows == 1;
    int d = mat->cols == vec->size;

    assert(a && b && !(c && d) || !(a && b) && c && d);

    if (a && b) {
        for (int i = 0; i < vec->size; i++)
            vec->array[i] = mat->array[i][0];
        return;
    }
    else if (c && d) {
        for (int i = 0; i < vec->size; i++)
            vec->array[i] = mat->array[0][i];
        return;
    }
}

int checkConvert(Matrix* mat, Vector* vec) {
    if (mat->rows == 1) {
        initVec(vec, mat->cols);
        for (int i = 0; i < mat->cols; i++)
            vec->array[i] = mat->array[0][i];
        return 1;
    }
    else if (mat->cols == 1) {
        initVec(vec, mat->rows);
        for (int i = 0; i < mat->rows; i++)
            vec->array[i] = mat->array[i][0];
        return 1;
    }
    
    return 0;
}

void convertVecToMat(Vector* vec, Matrix* mat) {
    int rows = vec->size == mat->rows;
    int cols = vec->size == mat->cols;
    int correct = mat->rows > 1 && mat->cols == 1 || mat->cols > 1 && mat->rows == 1;
    assert((rows || cols) && correct);

    if (rows) {
        for (int i = 0; i < vec->size; i++)
            mat->array[i][0] = vec->array[i];
    }
    else if (cols) {
        for (int i = 0; i < vec->size; i++)
            mat->array[0][i] = vec->array[i];
    }
}

//--------------------------------------------------------\\
// Arithmetic Operations                                  \\
//--------------------------------------------------------\\

void assignVecSparseRow(Vector* vec, int row, Sparse* sparse) {
    assert(sparse->cols == vec->size);

    for (int i = 0; i < sparse->arraySize; i++) {
        for (int j = 0; j < sparse->array[i].arraySize; j++) {
            if (sparse->array[i].array[j].first == row) {
                int c = sparse->array[i].value;
                int v = sparse->array[i].array[j].second;
                vec->array[c] = v;
            }
        }
    }
}

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

void scalarSparseMult(Sparse* mat, double scalar) {
    for (int i = 0; i < mat->arraySize; i++) {
        for (int j = 0; j < mat->array[i].arraySize; j++)
            mat->array[i].array[j].second *= scalar;
    }
}

void scalarSparseDiv(Sparse* mat, double scalar) {
    for (int i = 0; i < mat->arraySize; i++) {
        for (int j = 0; j < mat->array[i].arraySize; j++)
            mat->array[i].array[j].second /= scalar;
    }
}

void matAdd(Matrix* m1, Matrix* m2, Matrix* out) {
    assert(m1->rows == m2->rows);
    assert(m1->cols == m2->cols);
    assert(m1->rows == out->rows);
    assert(m2->cols == out->cols);

    for (int i = 0; i < m1->rows; i++) {
        for (int j = 0; j < m1->cols; j++)
            out->array[i][j] = m1->array[i][j] + m2->array[i][j];
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

// component function
// for use in matDot
// dot product of m1 row and m2 column
double matDotPartial(Matrix* m1, int row, Matrix* m2, int col) {
    double output = 0;
    for (int i = 0; i < m2->rows; i++)
        output += m1->array[row][i] * m2->array[i][col];

    return output;
}

double matVecDotPartialRow(Matrix* mat, int row, Vector* vec) {
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
        out->array[i] = matVecDotPartialRow(mat, i, vec);
}

void vecMatDot(Vector* vec, Matrix* mat, Vector* out) {
    assert(vec->size == mat->rows);
    assert(mat->cols == out->size);

    for (int i = 0; i < vec->size; i++)
        out->array[i] = matVecDotPartialCol(mat, i, vec);
}

// general matrix multiplication function
// not for use with vectors
void matDot(Matrix* m1, Matrix* m2, Matrix* m3) {
    Vector v1, v2, out;

    int m1Vec = m1->cols == 1 && m1->rows != 1 || m1->cols != 1 && m1->rows == 1;
    int m2Vec = m2->cols == 1 && m2->rows != 1 || m2->cols != 1 && m2->rows == 1;

    if (m1Vec && m2Vec) {
        assert(m3->rows == 1 && m3->cols == 1);

        initVec(&v1, m1->cols);
        initVec(&v2, m2->rows);
        convertMatToVec(m1, &v1);
        convertMatToVec(m2, &v2);
        
        m3->array[0][0] = vecDot(&v1, &v2);

        cleanVec(&v1);
        cleanVec(&v2);

        return;
    }

    else if (checkConvert(m1, &v1)) {
        initVec(&out, m2->cols);
        vecMatDot(&v1, m2, &out);
        convertVecToMat(&out, m3);
        cleanVec(&out);
        cleanVec(&v1);
        return;
    }

    else if (checkConvert(m2, &v1)) {
        initVec(&out, m1->rows);
        matVecDot(m1, &v1, &out);
        convertVecToMat(&out, m3);
        cleanVec(&out);
        cleanVec(&v1);
        return;
    }

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
void sparseVecDot(Sparse* sparse, Vector* vec, Vector* out) {
    assert(out->size == vec->size || out->size == sparse->cols || vec->size == sparse->rows);

    for (int i = 0; i < vec->size; i++) {
        double sum = 0;
        Vector row;
        initVec(&row, vec->size);
        assignVecSparseRow(&row, i, sparse);
        out->array[i] = vecDot(&row, vec);
        cleanVec(&row);
    }           
}

void vecSparseDot(Vector* vec, Sparse* sparse, Vector* out) {
    assert(out->size == sparse->cols);

    for (int i = 0; i < sparse->arraySize; i++) {
        double sum = 0;
        for (int j = 0; j < sparse->array[i].arraySize; j++) {
            double sparseValue = sparse->array[i].array[j].second;
            double vecValue = vec->array[sparse->array[i].array[j].first];
            
            sum += vecValue * sparseValue;
        }
        out->array[sparse->array[i].value] = sum;
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

void sparseDotFirst(Sparse* sparse, Matrix* mat, Matrix* out) {
    int matVecRow = mat->rows != 1 && mat->cols == 1;
    int matVecCol = mat->rows == 1 && mat->cols != 1;
    Vector vec, vOut;

    if (matVecRow) {
        initVec(&vec, mat->rows);
        initVec(&vOut, mat->rows);
        convertMatToVec(mat, &vec);
        sparseVecDot(sparse, &vec, &vOut);
        convertVecToMat(&vOut, out);
        cleanVec(&vec);
        cleanVec(&vOut);
        return;
    }
    else if (matVecCol) {
        initVec(&vec, mat->cols);
        initVec(&vOut, mat->cols);
        convertMatToVec(mat, &vec);
        sparseVecDot(sparse, &vec, &vOut);
        convertVecToMat(&vOut, out);
        cleanVec(&vec);
        cleanVec(&vOut);
        return;
    }

    assert(sparse->cols == mat->rows);
    assert(sparse->rows == out->rows);
    assert(mat->cols == out->cols);

    for (int i = 0; i < sparse->rows; i++) {
        Vector sRow;
        initVec(&sRow, sparse->cols);
        assignVecSparseRow(&sRow, i, sparse);
        for (int j = 0; j < mat->cols; j++)
            out->array[i][j] = matVecDotPartialCol(mat, j, &sRow);

        cleanVec(&sRow);
    }
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

// assign vector from matrix column
void assignVecMatCol(Vector* vec, int col, Matrix* mat) {
    assert(mat->rows == vec->size);

    for (int i = 0; i < mat->rows; i++)
        vec->array[i] = mat->array[i][col];
}

// append m2 to m1, row-wise
void appendMatRow(Matrix** m1, Matrix* m2) {
    assert((*m1)->cols == m2->cols);

    Matrix* pm1 = *m1;
    Matrix temp;
    initMat(&temp, pm1->rows + m2->rows, pm1->cols);

    for (int i = 0; i < pm1->rows; i++) {
        for (int j = 0; j < pm1->cols; j++)
            temp.array[i][j] = pm1->array[i][j];
    }

    for (int i = 0; i < m2->rows; i++) {
        for (int j = 0; j < m2->cols; j++)
            temp.array[i+pm1->rows][j] = m2->array[i][j];
    }

    cleanMat(pm1);
    *m1 = &temp;
}

void cloneVec(Vector* v1, Vector* v2) {
    assert(v1->size == v2->size);

    for (int i = 0; i < v1->size; i++)
        v2->array[i] = v1->array[i];
}

void cloneMat(Matrix* m1, Matrix* m2) {
    assert(m1->rows == m2->rows);
    assert(m1->cols == m2->cols);

    for (int i = 0; i < m1->rows; i++) {
        for (int j = 0; j < m1->cols; j++)
            m2->array[i][j] = m1->array[i][j];
    }
}

// reshape but only shrinking dimensions
void shrinkMat(Matrix* mat, int rows, int cols) {
    double** newArray = calloc(rows, sizeof(*newArray));
    for (int i = 0; i < rows; i++)
        newArray[i] = calloc(cols, sizeof(*newArray[i]));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++)
            newArray[i][j] = mat->array[i][j];
    }
    
    cleanMat(mat);
    
    mat->array = newArray;
    mat->rows = rows;
    mat->cols = cols;
}

void sparseToMat(Sparse* sparse, Matrix* mat) {
    assert(sparse->rows == mat->rows);
    assert(sparse->cols == mat->cols);
    
    for (int i = 0; i < sparse->arraySize; i++) {
        for (int j = 0; j < sparse->array[i].arraySize; j++)
            mat->array[sparse->array[i].array[j].first][sparse->array[i].value] = sparse->array[i].array[j].second;
    }
}

#endif
