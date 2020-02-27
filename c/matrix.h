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
    float* ptr;
    int rows;
    int cols;
} Matrix;

// TODO: Organizing this file
// TODO: Eigenvalues
// TODO: Formatting - sparseDotPartial looks like a disaster

//----------------------------------------\\
// Initialization                         \\
//----------------------------------------\\

void initMat(Matrix* mat, int r, int c) {
    mat->ptr = malloc(sizeof mat->ptr * r * c);
    mat->rows = r;
    mat->cols = c;
}

void clean(Matrix* mat) {
    free(mat->ptr);
}

void initRandom(Matrix* mat, int r, int c) {
    mat->ptr = malloc(sizeof mat->ptr * r * c);
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->ptr[i*mat->cols+j] = randomFloat();
    }
}

// uses normally distributed values
// uses ParentArray as a list of columns
// each nested Parent acting as a unique column
// with a list of rows in which there is a value
// for each row where there is a value in the sparse matrix
// eliminates unnecessary storage of zeros
void initSparse(ParentArray* mat, int r, int c, float density) {
    mat->array = malloc(sizeof mat->array);
    mat->arraySize = 1;
    mat->rows = r;
    mat->cols = c;
    
    int count = 0;
    int total = (int)(r*c*density);

    while (count < total) {
        int column = randRange(c);
        int in = parentArrayIn(mat, column);
        if (in+1 != mat->arraySize) {
            int row = randRange(r);
            int in2 = parentInFirst(&(mat->array[in]), row);
            if (in2 == 0) {
                parentAdd(&(mat->array[in]), row, marsagliaPolar());
                count++;
            }
        }
        else {
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
    mat->ptr = malloc(sizeof mat->ptr * r * c);
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->ptr[i*mat->cols+j] = marsagliaPolar();
    }
}

// shapes must be correct before use
void initTranspose(Matrix* mat, Matrix* matT) {
    assert(matT->rows != mat->cols || matT->cols != mat->rows);

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            matT->ptr[j*matT->cols+i] = mat->ptr[i*mat->cols+j];
    }
}

void printMat(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        printf("[ ");
        for (int j = 0; j < mat->cols; j++)
            printf("%f ", mat->ptr[i*mat->cols+j]);
        printf("]\n");
    }
}

void set(Matrix* mat, float value, int row, int col) {
    mat->ptr[row*mat->cols+col] = value;
}

//--------------------------------------------------------\\
// Arithmetic Operations                                  \\
//--------------------------------------------------------\\

float matDotPartial(Matrix* m1, int row, Matrix* m2, int col) {
    float output = 0;
    for (int i = 0; i < m2->rows; i++)
        output += m1->ptr[row*m1->cols+i] * m2->ptr[i*m2->cols+col];

    return output;
}

void matDot(Matrix* m1, Matrix* m2, Matrix* m3) {
    assert(m3->rows != m1->rows || m3->cols != m2->cols || m1->cols != m2->rows);

    for (int i = 0; i < m1->rows; i++) {
        for (int j = 0; j < m2->cols; j++)
            m3->ptr[i*m3->cols+j] = matDotPartial(m1, i, m2, j);
    }
}

float sparseDotPartial(Matrix* mat, int row, ParentArray* sparse, int col) {
    float output = 0;
    for (int i = 0; i < sparse->array[col].arraySize; i++)
        output += mat->ptr[row*mat->cols+sparse->array[col].array[i].first] * sparse->array[col].array[i].second;

    return output;
}

void sparseDotSecond(Matrix* mat, ParentArray* sparse, Matrix* out) {
    assert(out->rows == mat->rows || out->cols == sparse->cols || mat->cols == sparse->rows);
    
    int sparseIndex = 0;
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < sparse->arraySize; j++)
            out->ptr[i*out->cols+sparse->array[j].value] = sparseDotPartial(mat, i, sparse, j);
    }
}


//--------------------------------------------------------\\
// Manipulation                                           \\
//--------------------------------------------------------\\

// uses power iteration to find the max eigenvalue
//float maxEigen(Matrix

#endif
