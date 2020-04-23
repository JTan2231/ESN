#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "linalg.h"
#include "generation.h"

int main() {
    srand(time(0));

    Sparse sparse;
    initSparse(&sparse, 20, 20, 0.05);
    Matrix mat, H, Q;
    initSparseToMat(&sparse, &mat);

    initMat(&H, 0, 0);
    initMat(&Q, mat.rows, 0);

    Vector f;
    initVecRandom(&f, mat.rows);

    kStepArnoldi(&mat, &Q, &H, &f, 0, mat.rows);

    printMat(&H);

    qrHess(&H);

    printMat(&H);
}
