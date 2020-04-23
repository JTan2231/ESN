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

#define TOL 0.0000000000009

// !! TODO: clean your room !!
// TODO: don't forget assertions
// TODO: make Arnoldi sparse functions actually utilize sparse struct

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

//--------------------------------------------------------\\
// Eigenvectors and Eigenvalues                           \\
//--------------------------------------------------------\\

// performs QR factorization using the Modified Gram-Schmidt
// NOTE: mat is overwritten here
// see Matrix Computation, page 232 - Golub, Van Loan (1996)
int qrFact(Matrix* mat, Matrix* Q, Matrix* R) {
    assert(mat->rows == Q->rows);
    assert(mat->cols == Q->cols);
    assert(mat->cols == R->cols);
    assert(mat->rows == mat->cols);
    assert(R->rows == R->cols);

    for (int i = 0; i < mat->cols; i++) {
        Vector v;
        initVecCol(mat, i, &v);
        R->array[i][i] = magnitude(&v);
        for (int j = 0; j < Q->rows; j++)
            Q->array[j][i] = v.array[j] / R->array[i][i];
        cleanVec(&v);

        for (int j = i+1; j < mat->cols; j++) {
            Vector q, a;
            initVecCol(Q, i, &q);
            initVecCol(mat, j, &a);
            R->array[i][j] = vecDot(&q, &a);

            for (int k = 0; k < mat->rows; k++)
                mat->array[k][j] = mat->array[k][j] - q.array[k]*R->array[i][j];

            cleanVec(&q);
            cleanVec(&a);
        }
    }
}

// dimensions:
// -- mat = (n, n)
// -- Q = (n, k)
// -- H = (k, k)
// -- f = (n)
//
// extends an Arnoldi factorization from size k to size k + p
// k == 0, f == random vector for initial factorization
//
// reshapes:
// -- Q
// -- -- (n, k) to (n, k+p)
// -- H
// -- -- (k, k) to (k+p, k+p)
//
// see Algorithm 2.1, page 3 of:
// -- https://pdfs.semanticscholar.org/9b61/ec78bb1605143b60a7aafcde9fa39cb7e14a.pdf
void kStepArnoldi(Matrix* mat, Matrix* Q, Matrix* H, Vector* f, int k, int p) {
    assert(Q->rows == mat->rows);
    assert(Q->cols == k);
    assert(H->rows == H->cols);
    assert(H->rows == k);
    assert(f->size == mat->rows);

    for (int i = 0; i < p; i++) {
        double beta = magnitude(f);
        if (beta < TOL) {
            printf("ZERO MAGNITUDE AT STEP %d. RETURNING.\n", i);
            return;
        }

        Vector qj;
        initVec(&qj, f->size);
        scalarVecDivOut(f, beta, &qj);
        appendMatColVec(Q, &qj);

        Vector w;
        initVec(&w, mat->rows);
        matVecDot(mat, &qj, &w);
        
        Matrix qT;
        initTranspose(Q, &qT);
        shrinkMat(&qT, qT.rows-1, qT.cols);

        Vector hj;
        initVec(&hj, qT.rows);
        matVecDot(&qT, &w, &hj);

        double alpha = vecDot(&qj, &w);

        growMat(H, H->rows+1, H->cols+1);

        if (k == 0 && i == 0)
            H->array[0][0] = beta;
        else
            H->array[H->rows-1][H->cols-2] = beta; // newest subdiagonal

        for (int j = 0; j < hj.size; j++)
            H->array[j][H->cols-1] = hj.array[j];

        H->array[H->rows-1][H->cols-1] = alpha;

        Vector vh;
        initVec(&vh, f->size);
        Matrix qClone;
        initClone(Q, &qClone);
        shrinkMat(&qClone, qClone.rows, qClone.cols-1);
        matVecDot(&qClone, &hj, &vh);

        zeroVec(f);
        vecSubOut(&w, &vh, f);
        scalarVec(&qj, alpha);
        vecSub(f, &qj);

        cleanVec(&qj);
        cleanVec(&w);
        cleanMat(&qT);
        cleanVec(&hj);
        cleanVec(&vh);
        cleanMat(&qClone);
    }
}

// double-shift Hessenberg QR algorithm
// see Algorithm 4.5 of:
// https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter4.pdf
int qrHess(Matrix* H) {
    int p = H->rows-1;
    int q;
    while (p > 1) {
        q = p-1;
        double s = H->array[q][q] + H->array[p][p];
        double t = H->array[q][q]*H->array[p][p] - H->array[q][p]*H->array[p][q];

        // first three elements of column M
        double x = pow(H->array[0][0], 2) + H->array[0][1]*H->array[1][0] - s*H->array[0][0] + t;
        double y = H->array[1][0] * (H->array[0][0] + H->array[1][1] - s);
        double z = H->array[1][0] * H->array[2][1];
        if (y < TOL && z < TOL) {
            printf("y and z < TOL\n");
            break;
        }
        if (s < TOL || t < TOL) {
            printf("s or t < TOL\n");
            break;
        }

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
            if (isnan(d)) break;

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
        assert(!isnan(u.array[0]));
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

// finds the magnitude of the spectral radius of given Hessenberg matrix
double specHess(Matrix* H) {
    assert(!qrHess(H));

    double spec = 0;
    for (int i = 0; i < H->rows-1; i++) {
        double w = H->array[i][i];
        double x = H->array[i][i+1];
        double y = H->array[i+1][i];
        double z = H->array[i+1][i+1];

        double a = 1;
        double b = w + z;
        double c = w*z - y*x;

        double root = b*b - 4*a*c;
        double eigP = (b + sqrt(root))/(2*a);
        double eigN = (b - sqrt(root))/(2*a);
        if (root > 0) {
            if (fabs(eigP) > fabs(spec))
                spec = eigP;
            if (fabs(eigN) > fabs(spec))
                spec = eigN;
        }
    }

    return spec;
}

// TODO: make Arnoldi utilize sparse struct
// TODO: use IRA
double spectralRadius(Sparse* sparse) {
    Matrix mat, Q, H;
    initSparseToMat(sparse, &mat);
    initMat(&Q, mat.rows, 0);
    initMat(&H, 0, 0);
    Vector f;
    initVecRandom(&f, mat.rows);

    kStepArnoldi(&mat, &Q, &H, &f, 0, mat.rows);

    double spec = specHess(&H);

    cleanMat(&mat);
    cleanMat(&Q);
    cleanMat(&H);
    cleanVec(&f);

    return spec;
}

// Implicitly Restarted Arnoldi (IRA)
// see: http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter11.pdf
// TODO: make
// TODO: generalize (don't use only spectral radius)
void implicitArnoldiSparse(Sparse* sparse) {}

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
void inverse(Matrix* mat, Matrix* inverse) {
    assert(mat->rows == mat->cols);

    Matrix lower, upper, d, eye;
    initMat(&lower, mat->rows, mat->cols);
    initMat(&upper, mat->rows, mat->cols);
    initMat(&d, mat->rows, mat->cols);
    initIdent(&eye, mat->rows, mat->cols);

    crout(mat, &lower, &upper);

    for (int column = 0; column < d.cols; column++) {
        for (int i = 0; i < d.rows; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += lower.array[i][j] * d.array[j][column];

            d.array[i][column] = (eye.array[i][column]-sum) / lower.array[i][i];
        }
    }

    // solve for x (the inverse of mat)
    for (int column = 0; column < inverse->cols; column++) {
        for (int i = inverse->rows-1; i >= 0; i--) {
            double sum = 0;
            for (int j = i; j < inverse->rows; j++)
                sum += upper.array[i][j] * inverse->array[j][column];

            inverse->array[i][column] = d.array[i][column] - sum;
        }
    }
}

#endif
