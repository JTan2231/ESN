/*
 * Author: Joey Tan
 * Date Created: 3-7-20
 * Last Edit: 3-7-20, Joey Tan
 */

#ifndef LINALG
#define LINALG
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "map.h"

//--------------------------------------------------------\\
// Basic Vector Operations                                \\
//--------------------------------------------------------\\

// normalizes the given vector
void normalize(Vector* vec, double magnitude) {
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
    for (int i = 0; i < vec->size; i++) {
        printf("Current square value: %.16lf\n", vec->array[i]);
        output += vec->array[i] * vec->array[i];
    }

    printf("Sum of squares: %.16lf\n", output);

    return sqrt(output);
}

//--------------------------------------------------------\\
// Eigenvectors and Eigenvalues                           \\
//--------------------------------------------------------\\

// see:
//  - http://www.cs.cmu.edu/afs/cs/academic/class/15859n-f16/Handouts/TrefethenBau/ArnoldiIteration-33.pdf
//  - https://en.wikipedia.org/wiki/Arnoldi_iteration
// a dense version because I'm an idiot
void arnoldiDense(Matrix* mat, Matrix* Q, Matrix* H) {
    assert(mat->rows == Q->rows);
    assert(mat->cols == H->cols);
    assert(mat->rows == mat->cols);
 
    Vector b;
    Vector q;
    initVecRandom(&b, mat->cols);
    initVec(&q, b.size);
    normalizeOut(&b, &q, magnitude(&b));
    for (int i = 0; i < mat->rows; i++) {
        assignMatColVec(Q, i, &q);
        Vector temp;
        initVec(&temp, mat->rows);
        cloneVec(&q, &temp);
        matVecDot(mat, &temp, &q);

        for (int j = 0; j <= i; j++) {
            Vector qj;
            initVecCol(Q, j, &qj);
            H->array[j][i] = vecDot(&qj, &q);
            scalarVec(&qj, H->array[j][i]);
            vecSub(&q, &qj);
            cleanVec(&qj);
        }

        double mag = magnitude(&q);
        H->array[i][i-1] = mag;
        
        if (mag < 0.00000005) {
            printf("ZERO MAGNITUDE\n");
            return;
        }
        scalarVecDivOut(&q, mag, &q);

        cleanVec(&temp);
    }

    cleanVec(&b);
    cleanVec(&q);

}

void arnoldiSparse(Sparse* sparse, Matrix* Q, Matrix* H) {
    assert(sparse->rows == Q->rows);
    assert(sparse->cols == H->cols);
    assert(sparse->rows == sparse->cols);

    Vector b;
    Vector q;
    initVecRandom(&b, sparse->cols);
    initVec(&q, b.size);
    normalizeOut(&b, &q, magnitude(&b));
    assignMatColVec(Q, 0, &q);
    for (int i = 1; i < sparse->rows; i++) {
        cleanVec(&q);
        initVec(&q, b.size);
        Vector temp;
        initVec(&temp, sparse->rows);
        assignVecMatCol(&temp, i-1, Q);
        sparseVecDot(sparse, &temp, &q);

        for (int j = 0; j <= i; j++) {
            printf("i, j: %d, %d\n", i, j);
            Vector qj;
            initVecCol(Q, j, &qj);
            H->array[j][i] = vecDot(&qj, &q);
            printf("q:\n");
            printVec(&q);
            printf("qj:\n");
            printVec(&qj);
            printf("H[j][i]: %.16lf\n", H->array[j][i]);
            scalarVec(&qj, H->array[j][i]);
            vecSub(&q, &qj);
            cleanVec(&qj);
        }
        
        double mag = magnitude(&q);
        H->array[i][i-1] = mag;
        
        if (mag < 0.0000000012) {
            printf("ZERO MAGNITUDE\n");
            // if Q isn't square 
            if (i != sparse->rows - 1)
                shrinkMat(Q, sparse->rows, i+1);
            return;
        }
        
        scalarVecDivOut(&q, mag, &q);
        assignMatColVec(Q, i, &q);

        cleanVec(&temp);
    }

    cleanVec(&b);
    cleanVec(&q);
}

void qrSparse(Matrix* hess, Matrix* orth, Matrix* out) {
    assert(out->rows == hess->rows);// == orth->rows);
    assert(out->cols == orth->cols);// == orth->cols);
    
    int iterations = 1000;
    Matrix temp, orthT;
    
    initTranspose(orth, &orthT);
    initMat(&temp, orth->cols, hess->cols);
        
    matDot(&orthT, hess, &temp);
    matDot(&temp, orth, out);
    
    cleanMat(&temp);
    
    for (int i = 0; i < iterations-1; i++) {
        initMat(&temp, orth->rows, hess->cols);
        
        matDot(&orthT, out, &temp);
        matDot(&temp, orth, out);
        
        cleanMat(&temp);
    }
}

void hessenbergSparse(Matrix* arn, Sparse* sparse, Matrix* out) {
    assert(sparse->cols == arn->rows);
    
    Matrix temp;
    Matrix arnT;
    initMat(&temp, sparse->rows, sparse->cols);
    initTranspose(arn, &arnT);
    sparseDotSecond(&arnT, sparse, &temp);
    matDot(&temp, arn, out);
    cleanMat(&temp);
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

//--------------------------------------------------------\\
// Matrix Inverse                                         \\
//--------------------------------------------------------\\

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
