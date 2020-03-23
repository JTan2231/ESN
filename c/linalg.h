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
// a dense version because I'm an ~~idiot~~ overachiever
void arnoldiDense(Matrix* mat, Matrix* Q, Matrix* H) {
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
            assert(i == mat->rows-1);
            printf("ZERO MAGNITUDE\n");
            return;
        }
        
        scalarVecDivOut(&v, mag, &q);
        cleanVec(&v);
    }
}

void arnoldiSparse(Sparse* sparse, Matrix* Q, Matrix* H) {
    Matrix A;
    initMat(&A, sparse->rows, sparse->cols);
    sparseToMat(sparse, &A);
    arnoldiDense(&A, Q, H);
}

void qrHess(Matrix* H, Matrix* out) {
    int o = 0;
    int p = H->rows-1;
    int q;
    while (p > 1) {
        /*if (o == 1) {
            printf("After Second Iteration:\n");
            printMat(H);
            return;
        }*/
        printf("p: %d\n", p);
        q = p-1;
        double s = H->array[q][q] + H->array[p-1][p-1];
        double t = H->array[q][q]*H->array[p][p] - H->array[q][p]*H->array[p][q];

        // first three elements of column M
        double x = powf(H->array[0][0], 2) + H->array[0][1]*H->array[1][0] - s*H->array[0][0] + t;
        double y = H->array[1][0] * (H->array[0][0] + H->array[1][1] - s);
        double z = H->array[1][0] * H->array[2][1];

        Vector u;
        Matrix temp;

        for (int i = -1; i < p-2; i++) {
            
            // determine Householder reflector 
            // represented only by vector u
            initVec(&u, 3);
            int rho = -1*sgn(x);
            double uMag = sqrt(x*x + y*y + z*z);
            double xTemp = x - rho*uMag;
            double d = sqrt(xTemp*xTemp + y*y + z*z);

            if (d > 1000) {
                printf("\nDIVERGING\n");
                return;
            }

            u.array[0] = (x - rho*uMag) / d;
            u.array[1] = y / d;
            u.array[2] = z / d;

            int r = max(0, i);

            /*printf("u:\n");
            printVec(&u);
            printf("u check:\n");
            if (1) {
                double ux = 2*(u.array[0]*x + u.array[1]*y + u.array[2]*z);
                printf("%lf\n%lf\n%lf\n", x-u.array[0]*ux, y-u.array[1]*ux, z-u.array[2]*ux);
            }

            printf("d: %lf\n", d);*/

            initMat(&temp, 3, H->rows-r+1);
            //for (int l = 0; l < 3; l++) {
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

                    /*if (l == 0) {
                        a -= 1 - u.array[0]*du;
                        b *= 0 - u.array[1]*du;
                        c *= 0 - u.array[2]*du;
                    }
                    else if (l == 1) {
                        a *= 0 - u.array*du;
                        b *= 1 - 2*du;
                        c *= 0 - 2*du;
                    }
                    else if (l == 2) {
                        a *= 0 - 2*du;
                        b *= 0 - 2*du;
                        c *= 1 - 2*du;
                    }

                    temp.array[l][k-r] = a + b + c;*/
                    //printf("First temp assignment: %lf, %lf, %lf\n", a-aTemp, b-bTemp, c-cTemp);
                }
            //}

            // move from temp to H
            /*for (int j = r; j < H->cols; j++) {
                H->array[i+1][j] = temp.array[0][j-r];
                H->array[i+2][j] = temp.array[1][j-r];
                H->array[i+3][j] = temp.array[2][j-r];
            }*/

            cleanMat(&temp);

            //printf("Right Householder:\n");
            //printMat(H);

            r = min(i+5, p+1);
            //printf("r: %d\n", r);

            initMat(&temp, r, 3);
            for (int l = 0; l < r; l++) {
                //for (int k = 0; k < 3; k++) {
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

                    /*if (k == 0) {
                        a *= 1 - 2*du;
                        b *= 0 - 2*du;
                        c *= 0 - 2*du;
                    }
                    else if (k == 1) {
                        a *= 0 - 2*du;
                        b *= 1 - 2*du;
                        c *= 0 - 2*du;
                    }
                    else if (k == 2) {
                        a *= 0 - 2*du;
                        b *= 0 - 2*du;
                        c *= 1 - 2*du;
                    }
                    
                    temp.array[l][k] = a + b + c;*/
                    //printf("Second temp assignment: %lf, %lf, %lf\n", a-aTemp, b-bTemp, c-cTemp);
                //}
            }

            // move from temp to H
            /*for (int j = 0; j < temp.rows; j++) {
                H->array[j][i+1] = temp.array[j][0];
                H->array[j][i+2] = temp.array[j][1];
                H->array[j][i+3] = temp.array[j][2];
            }*/

            cleanMat(&temp);
            cleanVec(&u);

            x = H->array[i+2][i+1];
            y = H->array[i+3][i+1];
            if (i < p-3)
                z = H->array[i+4][i+1];

            //printf("Left Householder:\n");
            //printMat(H);
            //return;
            /*if (i == 0) {
                printf("Two Householder cycles:\n");
                printMat(H);
                return;
            }*/
        }

        initVec(&u, 2);
        double uMag = sqrt(x*x + y*y);
        int rho = -1*sgn(x);
        double xTemp = x - rho*uMag;
        double d = sqrt(xTemp*xTemp + y*y);
        u.array[0] = (x - rho*uMag) / d;
        u.array[1] = y / d;

        /*printf("Givens u check:\n");
        if (1) {
            double ux = 2*(u.array[0]*x + u.array[1]*y);
            printf("%lf\n%lf\n", x-u.array[0]*ux, y-u.array[1]*ux);
        }*/

        // Givens rotation P from the left
        initMat(&temp, 2, H->cols-p+2);
        //for (int j = 0; j < 2; j++) {
            for (int k = p-2; k < H->cols; k++) {
                double a = H->array[q][k];
                double b = H->array[p][k];

                double ux = 2*(u.array[0]*a + u.array[1]*b);

                double aTemp = u.array[0] * ux;
                double bTemp = u.array[1] * ux;

                H->array[q][k] = a - aTemp;
                H->array[p][k] = b - bTemp;

                /*if (j == 0) {
                    a *= 1 - 2*du;
                    b *= 0 - 2*du;
                }
                else if (j == 1) {
                    a *= 0 - 2*du;
                    b *= 1 - 2*du;
                }
                
                temp.array[j][k-p+2] = a + b;*/
                //printf("Third temp assignment: %lf, %lf\n", a-aTemp, b-bTemp);
            }
        //}

        // move temp to H
        /*for (int j = q; j < p+1; j++) {
            for (int k = p-2; k < H->cols; k++)
                H->array[j][k] = temp.array[j-q][k-p+2];
        }*/

        cleanMat(&temp);

        initMat(&temp, p, 2);
        for (int j = 0 ; j < p; j++) {
            //for (int k = p-1; k < p+1; k++) {
                double a = H->array[j][p-1];
                double b = H->array[j][p];

                double ux = 2*(u.array[0]*a + u.array[1]*b);

                double aTemp = u.array[0] * ux;
                double bTemp = u.array[1] * ux;

                H->array[j][p-1] = a - aTemp;
                H->array[j][p] = b - bTemp;

                //printf("Fourth temp assignment: %lf, %lf\n", a-aTemp, b-bTemp);
            //}
        }

        // move temp to H
        /*for (int j = 0; j < p; j++) {
            for (int k = p-1; k < p+1; k++)
                H->array[j][k] = temp.array[j][k-p+1];
        }*/
        
        if (isnan(H->array[p][q])) {
            printf("NAN\n");
            return;
        }
        printf("H->array[%d][%d]: %.16lf\n", p, q, H->array[p][q]);
        if (fabs(H->array[p][q]) < TOL*(fabs(H->array[q][q]) + fabs(H->array[p][p]))) {
            printf("Zero check 1:\n");
            H->array[p][q] = 0;
            p--;
            q = p - 1;
        }
        else if (fabs(H->array[p-1][q-1]) < TOL*(fabs(H->array[q-1][q-1]) + fabs(H->array[q][q]))) {
            printf("Zero check 2:\n");
            H->array[p-1][q-1] = 0;
            p -= 2;
            q = p - 1;
        }

        //printf("Post-Householders check:\n");
        //printMat(H);
        //return;

        o++;
    }

    printf("Final:\n");
    printMat(H);
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
