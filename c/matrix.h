/*
 * Author: Joey Tan
 * Date Created: 2-18-20
 * Last Edit: 2-22-20, Joey Tan
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "generation.h"

typedef struct {
    float* ptr;
    int rows;
    int cols;
} Matrix;

// TODO: Organizating this file

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
//void initSparse(Matrix* mat, int r, int c, float density) {


// initializes a matrix using normally distributed random values
void initRandomNormal(Matrix* mat, int r, int c) {
    mat->ptr = malloc(sizeof mat->ptr * r * c);
    mat->rows = r;
    mat->cols = c;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            mat->ptr[i*mat->cols+j] = marsagliaNormal();
    }
}

// shapes must be correct before use
void initTranspose(Matrix* mat, Matrix* matT) {
    if (matT->rows != mat->cols || matT->cols != mat->rows) {
        printf("ERROR: Matrix shapes incompatible. Transpose incomplete.\n");
        return;
    }

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
    if (m3->rows != m1->rows || m3->cols != m2->cols || m1->cols != m2->rows) {
        printf("ERROR: Matrix shapes incompatible. Required shapes (rows, columns): (M, N) x (N, P) -> (M, P)\n");
        return;
    }

    for (int i = 0; i < m1->rows; i++) {
        for (int j = 0; j < m2->cols; j++)
            m3->ptr[i*m3->cols+j] = matDotPartial(m1, i, m2, j);
    }
}

