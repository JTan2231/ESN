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
#include <matrix.h>
#include <sparse.h>

// !! TODO: clean your room !!
// TODO: don't forget assertions
// TODO: make Arnoldi sparse functions actually utilize sparse struct

//--------------------------------------------------------\\
// Basic Vector Operations                                \\
//--------------------------------------------------------\\

void normalize(Vector* vec, double magnitude);
void normalizeOut(Vector* vec, Vector* output, double magnitude);
double magnitude(Vector* vec);

//--------------------------------------------------------\\
// Miscellaneous Math Operations                          \\
//--------------------------------------------------------\\

// returns the sign (-1, 0, 1) of x
int sgn(double x);
int max(int a, int b);
int min(int a, int b);

//--------------------------------------------------------\\
// Eigenvalues                                            \\
//--------------------------------------------------------\\

// computes QR Factorization of matrix mat, with output in matrices Q and R
int qrHess(Matrix* H);

// expands Arnoldi factorization from size k to size k+p
// if k == 0 (i.e. initial factorization), f should be a random vector
void kStepArnoldi(Matrix* mat, Matrix* Q, Matrix* H, Vector* f, int k, int p);
double specHess(Matrix* H);
double spectralRadius(Sparse* sparse);

//--------------------------------------------------------\\
// Matrix Inverse                                         \\
//--------------------------------------------------------\\

void crout(Matrix* mat, Matrix* lower, Matrix* upper);
void inverse(Matrix* mat, Matrix* inverse);

#endif
