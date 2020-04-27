#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "linalg.h"
#include "generation.h"

int main() {
    srand(time(0));

    Sparse sparse;
    initSparse(&sparse, 6, 6, 0.2);
    Matrix mat, H, Q;
    initSparseToMat(&sparse, &mat);
    initMat(&Q, sparse.rows, 0);
    initMat(&H, 0, 0);
    Vector f;
    initVecRandom(&f, sparse.rows);

    kStepArnoldi(&mat, &Q, &H, &f, 0, sparse.rows);

    qrHess(&H);

    eigsHess(&H);

    printf("Original:\n");
    printMat(&mat);

    //double spec = spectralRadius(&sparse);

    //printMat(&H);
    //printf("spectral radius: %.15lf\n", spec);
}
