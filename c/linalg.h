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

// TODO: clean your room
// TODO: don't forget assertions
// TODO: make Arnoldi sparse functions actually utilize sparse struct

//--------------------------------------------------------\\
// Basic Vector Operations                                \\
//--------------------------------------------------------\\

Complex compRealDiv(Complex a, double b) {
    Complex out;
    out.real = a.real / b;
    out.i = a.i / b;

    return out;
}

Complex compAdd(Complex a, Complex b) {
    Complex out;
    out.real = a.real + b.real;
    out.i = a.i + b.i;

    return out;
}

Complex compSub(Complex a, Complex b) {
    Complex out;
    out.real = a.real - b.real;
    out.i = a.i - b.i;

    return out;
}

Complex compMult(Complex a, Complex b) {
    Complex out;
    out.real = a.real * b.real;
    out.i = a.i - b.i;

    return out;
}

Complex compScalarMult(Complex a, double b) {
    Complex out;
    out.real = a.real * b;
    out.i = a.i * b;

    return out;
}

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

// calculates the magnitude of given complex vector
double complexMagnitude(ComplexVector* compVec) {
    double output = 0;
    for (int i = 0; i < compVec->size; i++) {
        double real = compVec->array[i].real;
        double comp = compVec->array[i].i;
        output += sqrt(real*real + comp*comp);
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
        if (mag < TOL) {
            if (i == mat->rows-1)
                return 0;

            shrinkMat(H, i, i);
            return 1;
        }

        scalarVecDivOut(&v, mag, &q);
        cleanVec(&v);
    }

    return 0;
}

int iStepArnoldiDense(Matrix* mat, Matrix* Q, Matrix* H, Vector* residual, int steps, int useRand, int start) {
    printf("iSTEP ARNOLDI DENSE:\n");
    assert(mat->rows == Q->rows);
    /*assert(Q->cols == steps);
    assert(H->rows == steps+1);
    assert(H->cols == steps);*/
    assert(residual->size = mat->rows);

    printf("-- residual:\n");
    printVec(residual);

    Vector q, v;
    if (useRand) {
        Vector b;
        initVec(&q, mat->rows);
        initVecRandomNormal(&b, mat->rows);
        normalizeOut(&b, &q, magnitude(&b));
    }
    else initVecClone(residual, &q);

    for (int i = start; i < steps; i++) {
        printf("-- step %d\n", i);
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
        if (mag < TOL) {
            printf("ZERO MAGNITUDE iSTEP\n");
            if (i == steps-1)
                return 0;

            shrinkMat(H, i, i);
            printf("ARNOLDI WARNING: H shrunk\n");
            return 1;
        }

        scalarVecDivOut(&v, mag, &q);
        cloneVec(&v, residual);
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

int iStepArnoldiSparse(Sparse* sparse, Matrix* Q, Matrix* H, Vector* residual, int steps, int useRand, int start) {
    Matrix A;
    initSparseToMat(sparse, &A);
    return iStepArnoldiDense(&A, Q, H, residual, steps, useRand, start);
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
        double x = pow(H->array[0][0], 2) + H->array[0][1]*H->array[1][0] - s*H->array[0][0] + t;
        double y = H->array[1][0] * (H->array[0][0] + H->array[1][1] - s);
        double z = H->array[1][0] * H->array[2][1];
        if (y == 0 && z == 0) break;

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
        printf("xTemp, y: %lf, %lf\n", xTemp, y);
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

// TODO: please format this monster
/*void qrHessComplex(ComplexMatrix* H) {
    int p = H->rows-1;
    int q;
    while (p > 1) {
        //printf("p: %d\n", p);
        q = p-1;
        Complex s = compAdd(H->array[q][q], H->array[p][p]);
        Complex t = compSub(compMult(H->array[q][q], H->array[p][p]), compMult(H->array[q][p], H->array[p][q]));

        // first three elements of column M
        Complex x = compSub(compAdd(compMult(H->array[0][0], H->array[0][0]), compMult(H->array[0][1], H->array[1][0])), compAdd(compMult(s*H->array[0][0]), t));
        Complex y = compAdd(compMult(H->array[1][0], (H->array[0][0])), compSub(H->array[1][1], s));
        Complex z = compMult(H->array[1][0], H->array[2][1]);

        ComplexVector u;

        for (int i = -1; i < p-2; i++) {    
            //printMat(H);
            // determine Householder reflector 
            // represented only by vector u
            initComplexVec(&u, 3);
            int rho = -1*sgn(x.real);
            double uMag = sqrt(sqrt(x.real*x.real + x.i*x.i) + sqrt(y.real*y.real + y.i*y.i) + sqrt(z.real*z.real + z.i*z.i));
            Complex xTemp = {xTemp.real -= rho*uMag, xTemp.i = x.i};
            double d = sqrt(sqrt(xTemp.real*xTemp.real + xTemp.i*xTemp.i) + sqrt(y.real*y.real + y.i*y.i) + sqrt(z.real*z.real + z.i*z.i));

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
}*/

int in(double* array, double value, int size) {
    for (int i = 0; i < size; i++) {
        if (array[i] == value)
            return 1;
    }

    return 0;
}

// returns index of value in array
// if no value then -1
int ind(double* array, double value, int size) {
    for (int i = 0; i < size; i++) {
        if (array[i] == value)
            return i;
    }

    return -1;
}

// finds the magnitude of the spectral radius of given Hessenberg matrix
double specHess(Matrix* H) {
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
            printf("CANDIDATES: %lf, %lf\n", eigP, eigN);
            if (fabs(eigP) > fabs(spec))
                spec = eigP;
            if (fabs(eigN) > fabs(spec))
                spec = eigN;
        }
    }

    return spec;
}

// gets a list of eigenvalues for a Hessenberg matrix
// ignores complex eigenvalues
double* eigsHess(Matrix* H) {
    Matrix workCopy;
    initClone(H, &workCopy);

    printf("EIGS HESS:\n");
    assert(!qrHess(&workCopy));
    printf("-- qrHess completed.\n");

    // solve for eigenvalues on the diagonal
    double* spec = calloc((workCopy.rows+1)*2, sizeof(double));
    spec[0] = workCopy.rows;
    int reals = 0;
    for (int i = 0; i < workCopy.rows-1; i++) {
        double w = workCopy.array[i][i];
        double x = workCopy.array[i][i+1];
        double y = workCopy.array[i+1][i];
        double z = workCopy.array[i+1][i+1];

        double a = 1;
        double b = w + z;
        double c = w*z - y*x;

        double root = b*b - 4*a*c;

        double solution = (b + sqrt(root)) / (2*a);

        if (!in(spec, solution, workCopy.rows)) {
            if (root < 0) {
                spec[i+1] = nan("");
                continue;
            }
            else {
                spec[i+1] = solution;
                reals++;
            }
        }
    }

    if (reals != workCopy.rows) {
        double* temp = calloc((reals+1)*2, sizeof(double));
        temp[0] = spec[0];
        int r = 1;
        double largest = 0;
        for (int i = 1; i < workCopy.rows+1; i++) {
            if (!isnan(spec[i])) {
                temp[r] = spec[i];
                r++;
            }
        }

        assert(r);

        free(spec);
        spec = temp;
    }

    printf("EIGS HESS FINISHED\n");
    return spec;
}

// Implicitly Restarted Arnoldi (IRA)
// see: http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter11.pdf
// TODO: generalize (don't use only spectral radius
void implicitArnoldiSparse(Sparse* sparse) {
    assert(sparse->rows == sparse->cols);

    int steps;
    if (sparse->rows < 6 || sparse->cols < 6)
        steps = sparse->rows - 2;
    else steps = 6;

    Matrix Q, H; // from m-step Arnoldi
    initMat(&Q, sparse->rows, steps);
    initMat(&H, steps+1, steps);

    Vector f; // residual
    initVec(&f, sparse->cols);

    if (iStepArnoldiSparse(sparse, &Q, &H, &f, steps, 1, 0))
        steps = H.rows;

    printf("iStep completed.\n");

    // TODO: convergence criterion
    for (int i = 0; i < 5; i++) {
        double* spec = eigsHess(&H);
        int range = (int)spec[0];
        double largest = 0;
        int specIndex = -1;
        for (int j = 1; j < range+1; j++) {
            if (fabs(spec[j] > largest)) {
                largest = fabs(spec[i]);
                specIndex = i;
            }
        }

        printf("IDENTITY INITIALIZATION\n");
        Matrix qTemp; // Q
        initIdent(&qTemp, steps, steps);
        printf("IDENTITY INITIALIZED\n");
        // implicit QR step
        for (int j = 1; j < range+1 && j != specIndex; j++) {
            Matrix shiftedH;
            initClone(&H, &shiftedH);
            shrinkMat(&shiftedH, H.rows-1, H.cols);
            for (int i = 0; i < shiftedH.rows; i++)
                shiftedH.array[i][i] -= spec[j];

            Matrix qFact, rFact; // QR factorization -- Q_j, R_j
            initMat(&qFact, shiftedH.rows, shiftedH.rows);
            initMat(&rFact, shiftedH.rows, shiftedH.cols);

            qrFact(&shiftedH, &qFact, &rFact);
            cleanMat(&rFact);
            cleanMat(&shiftedH);
            printf("-- QR factorization completed\n");

            shrinkMat(&H, H.rows-1, H.cols);

            // H_m = Q_j^T * H_m * Q_j
            Matrix temp, qFactT;
            initMat(&temp, qFact.rows, H.cols);
            initTranspose(&qFact, &qFactT);
            matDot(&qFactT, &H, &temp);
            zeroMat(&H);
            matDot(&temp, &qFact, &H);
            cleanMat(&temp);
            cleanMat(&qFactT);

            growMat(&H, H.rows+1, H.cols);

            // Q = Q * Q_j
            initMat(&temp, qTemp.rows, qFact.cols);
            matDot(&qTemp, &qFact, &temp);
            cloneMat(&temp, &qTemp);
            cleanMat(&temp);
        }

        double beta = H.array[1][0];
        double sigma = Q.array[steps-1][0];
        Vector v;
        initVecCol(&Q, 1, &v);
        scalarVec(&v, beta);
        scalarVec(&f, sigma);
        vecAdd(&f, &v);

        // "shrink" matrix qTemp
        for (int i = 0; i < qTemp.rows; i++) {
            for (int j = 1; j < qTemp.cols; j++)
                qTemp.array[i][j] = 0;
        }

        Matrix temp;
        initMat(&temp, Q.rows, Q.cols);
        matDot(&Q, &qTemp, &temp);
        cloneMat(&temp, &Q);
        cleanMat(&temp);

        // "shrink" matrix H
        for (int i = 1; i < H.rows; i++) {
            for (int j = 0; j < H.cols; j++)
                H.array[i][j] = 0;
        }

        for (int i = 1; i < H.cols; i++)
            H.array[0][i] = 0;

        iStepArnoldiSparse(sparse, &Q, &H, &f, steps, 0, 1);

        free(spec);
    }

    shrinkMat(&H, H.rows-1, H.cols);

    printf("FINISHED. H:\n");
    printMat(&H);
    double out = specHess(&H);
    printf("spectral radius: %lf\n", out);
}

/*void eigsHess(ComplexMatrix* H, Vector* eigs) {
    assert(eigs->array == NULL); // must be uninitialized

    ComplexMatrix workCopy;
    initComplexFromReal(H, &workCopy);

    qrHess(&workCopy);

    // solve for eigenvalues on the diagonal
    for (int i = 0; i < workCopy.rows-1; i++) {
        double w = workCopy.array[i][i];
        double x = workCopy.array[i][i+1];
        double y = workCopy.array[i+1][i];
        double z = workCopy.array[i+1][i+1];

        double a = 1;
        double b = w + z;
        double c = w*z - y*x;

        double root = b*b - 4*a*c;

*/

double spectralRadius(Sparse* sparse) {
    int STEPS = sparse->rows;
    Matrix Q, H, mat;
    initMat(&Q, sparse->rows, STEPS);
    initMat(&H, STEPS+1, STEPS);
    initMat(&mat, sparse->rows, sparse->cols);

    int o = 0;
    int threshold = 1000;
    /*while (arnoldiSparse(sparse, &Q, &H) != 0 && o < threshold) {
        reinitSparse(sparse, sparse->rows, sparse->cols, sparse->density);
        reinitMat(&Q, Q.rows, Q.cols);
        reinitMat(&H, H.rows, H.cols);
        o++;
    }*/

    arnoldiSparse(sparse, &Q, &H);

    if (o >= threshold) {
        printf("-- ARNOLDI FAILED. ABORTING.\n");
        assert(0);
    }

    printf("-- Arnoldi completed.\n");

    assert(!qrHess(&H));
    printf("-- QR Algorithm complete.\n");
    double spec = 0;
    for (int i = 0; i < H.rows; i++) {
        if (fabs(H.array[i][i]) > spec)
            spec = fabs(H.array[i][i]);
    }

    printf("H:\n");
    printMat(&H);

    return spec;

    printf("H:\n");
    printMat(&H);

    return 0;

    cleanMat(&Q);
    cleanMat(&H);
    cleanMat(&mat);
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
