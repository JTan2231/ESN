/*
 * Author: Joey Tan
 * Date Created: 3-7-20
 * Last Edit: 3-13-20, Joey Tan
 */

#ifndef LINALG
#define LINALG
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "map.h"

#define TOL 0.0000012

// TODO: clean your room
// TODO: don't forget assertions

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
        output += vec->array[i] * vec->array[i];
    }

    return sqrt(output);
}

double sign(double a, double b) {
    if (b < 0) return -1*a;
    return a;
}

int sgn(double x) {
    if (x < 0) return -1;
    if (x > 0) return 1;
    return 0;
}

int max(int a, int b) {
    if (a > b) return a;
    return b;
}

int min(int a, int b) {
    if (a < b) return a;
    return b;
}

void genGivens(double* ap, double* bp, double* cp, double* sp) {
    if (*bp != 0) {
        double r = hypot(*ap, *bp);
        *cp = *ap/r;
        *sp = -1*(*bp) * 1./r;
        return;
    }

    *cp = 1;
    *sp = 0;
}

// performs a Givens rotation on the given matrix (one iteration)
// note: this and givensRight are specifically for the QR algorithm
// this function builds rArray for use in givensRight()
void givensLeft(Matrix* mat, int iRow) {
    assert(mat->rows == mat->cols);

    double a, b, c, s;

    a = mat->array[iRow][iRow];
    b = mat->array[iRow+1][iRow];

    double *ap = &a, *bp = &b, *cp = &c, *sp = &s;

    genGivens(ap, bp, cp, sp);
    
    for (int i = iRow; i < mat->cols-1; i++) {
        a = mat->array[iRow][i];
        b = mat->array[iRow+1][i];

        double temp = c*a - s*b;
        mat->array[iRow+1][i] = c*b + s*a;
        mat->array[iRow][i] = temp;
    }
}

void givensRight(Matrix* mat, int iCol) {
    assert(mat->rows == mat->cols);
    
    double a, b, c, s;
    a = mat->array[0][iCol];
    b = mat->array[0][iCol+1];
    
    double *ap = &a, *bp = &b, *cp = &c, *sp = &s;

    genGivens(ap, bp, cp, sp);

    for (int i = 0; i < iCol+2; i++) {
        a = mat->array[i][iCol];
        b = mat->array[i][iCol+1];

        double temp = c*a - s*b;
        mat->array[i][iCol+1] = c*b + s*a;
        mat->array[i][iCol] = temp;
    }
}

//--------------------------------------------------------\\
// Eigenvectors and Eigenvalues                           \\
//--------------------------------------------------------\\

// see:
//  - http://www.cs.cmu.edu/afs/cs/academic/class/15859n-f16/Handouts/TrefethenBau/ArnoldiIteration-33.pdf
//  - https://en.wikipedia.org/wiki/Arnoldi_iteration
// a dense version because I like to characterize myself in my code
int arnoldiDense(Matrix* mat, Matrix* Q, Matrix* H) {
    assert(mat->rows == Q->rows);
    assert(mat->cols == H->rows-1);
    assert(mat->rows == mat->cols);

    Vector b, q, v;
    initVec(&q, mat->rows);
    initVecRandomNormal(&b, mat->rows);
    normalizeOut(&b, &q, magnitude(&b));
    for (int i = 0; i < mat->rows; i++) {
        initVec(&v, mat->rows);
        matVecDot(mat, &q, &v);
        assignMatColVec(Q, i, &q);
        
        for (int j = 0; j <= i; j++) {
            Vector qj;
            initVec(&qj, Q->rows);
            assignVecMatCol(&qj, j, Q);
            H->array[j][i] = vecDot(&qj, &v);
            scalarVec(&qj, H->array[j][i]);
            vecSub(&v, &qj);
            cleanVec(&qj);
        }
        
        double mag = magnitude(&v);
        H->array[i+1][i] = mag;
        if (mag < 0.00000000012) {
            //printf("i: %d\n", i);
            if (i == mat->rows-1)
                return 0;

            //printf("ZERO MAGNITUDE\n");
            return 1;
        }
        
        scalarVecDivOut(&v, mag, &q);
        cleanVec(&v);
    }

    return 0;
}

int arnoldiSparse(Sparse* sparse, Matrix* Q, Matrix* H) {
    Matrix A;
    initMat(&A, sparse->rows, sparse->cols);
    sparseToMat(sparse, &A);
    arnoldiDense(&A, Q, H);
}

int qrHess(Matrix* H) {
    int p = H->rows-1;
    int q;
    while (p > 1) {
        //printf("p: %d\n", p);
        q = p-1;
        double s = H->array[q][q] + H->array[p][p];
        double t = H->array[q][q]*H->array[p][p] - H->array[q][p]*H->array[p][q];

        // first three elements of column M
        double x = powf(H->array[0][0], 2) + H->array[0][1]*H->array[1][0] - s*H->array[0][0] + t;
        double y = H->array[1][0] * (H->array[0][0] + H->array[1][1] - s);
        double z = H->array[1][0] * H->array[2][1];

        Vector u;

        for (int i = -1; i < p-2; i++) {    
            //printMat(H);
            // determine Householder reflector 
            // represented only by vector u
            initVec(&u, 3);
            int rho = -1*sgn(x);
            double uMag = sqrt(x*x + y*y + z*z);
            double xTemp = x - rho*uMag;
            double d = sqrt(xTemp*xTemp + y*y + z*z);

            if (d > 1000) {
                printf("-- QR FAILED: DIVERGING\n");
                return 1;
            }

            u.array[0] = (x - rho*uMag) / d;
            u.array[1] = y / d;
            u.array[2] = z / d;

            int r = max(0, i);
            for (int k = r; k < H->cols; k++) {
                double a = H->array[i+1][k];
                double b = H->array[i+2][k];
                double c = H->array[i+3][k];

                double ux = 2*(u.array[0]*a + u.array[1]*b + u.array[2]*c);
                
                double aTemp = u.array[0] * ux;
                double bTemp = u.array[1] * ux;
                double cTemp = u.array[2] * ux;

                H->array[i+1][k] = a - aTemp;
                H->array[i+2][k] = b - bTemp;
                H->array[i+3][k] = c - cTemp;
            }

            r = min(i+4, p);
            for (int l = 0; l < r+1; l++) {
                double a = H->array[l][i+1];
                double b = H->array[l][i+2];
                double c = H->array[l][i+3];

                double ux = 2*(u.array[0]*a + u.array[1]*b + u.array[2]*c);

                double aTemp = u.array[0] * ux;
                double bTemp = u.array[1] * ux;
                double cTemp = u.array[2] * ux;

                H->array[l][i+1] = a - aTemp;
                H->array[l][i+2] = b - bTemp;
                H->array[l][i+3] = c - cTemp;
            }
            
            cleanVec(&u);

            x = H->array[i+2][i+1];
            y = H->array[i+3][i+1];
            if (i < p-3)
                z = H->array[i+4][i+1];
        }

        initVec(&u, 2);
        double uMag = sqrt(x*x + y*y);
        int rho = -1*sgn(x);
        double xTemp = x - rho*uMag;
        double d = sqrt(xTemp*xTemp + y*y);
        u.array[0] = (x - rho*uMag) / d;
        u.array[1] = y / d;

        // Givens rotation P from the left
        for (int k = p-2; k < H->cols; k++) {
            double a = H->array[q][k];
            double b = H->array[p][k];

            double ux = 2*(u.array[0]*a + u.array[1]*b);

            double aTemp = u.array[0] * ux;
            double bTemp = u.array[1] * ux;

            H->array[q][k] = a - aTemp;
            H->array[p][k] = b - bTemp;
        }

        for (int j = 0 ; j < p+1; j++) {
            double a = H->array[j][p-1];
            double b = H->array[j][p];

            double ux = 2*(u.array[0]*a + u.array[1]*b);

            double aTemp = u.array[0] * ux;
            double bTemp = u.array[1] * ux;

            H->array[j][p-1] = a - aTemp;
            H->array[j][p] = b - bTemp;
        }
        
        if (isnan(H->array[p][q])) {
            printf("-- QR FAILED: NAN\n");
            return 1;
        }
        
        if (fabs(H->array[p][q]) < TOL*(fabs(H->array[q][q]) + fabs(H->array[p][p]))) {
            H->array[p][q] = 0;
            p--;
            q = p - 1;
        }
        else if (fabs(H->array[p-1][q-1]) < TOL*(fabs(H->array[q-1][q-1]) + fabs(H->array[q][q]))) {
            H->array[p-1][q-1] = 0;
            p -= 2;
            q = p - 1;
        }
    }

    return 0;
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
