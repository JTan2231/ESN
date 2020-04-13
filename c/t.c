#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "linalg.h"
#include "generation.h"

int main() {
    srand(time(0));

    /*Matrix mat, Q, R;
    initRandom(&mat, 5, 5);
    initIdent(&Q, 5, 5);
    initMat(&R, 5, 5);

    qrFact(&mat, &Q, &R);

    printf("mat:\n");
    printMat(&mat);
    printf("Q:\n");
    printMat(&Q);
    printf("R:\n");
    printMat(&R);

    printf("Orthogonality check\n");

    Matrix qt, out;
    initTranspose(&Q, &qt);
    initMat(&out, 5, 5);
    matDot(&Q, &qt, &out);

    printMat(&out);*/

    Sparse sparse;
    initSparse(&sparse, 6, 6, 0.4);
    Matrix mat;
    initSparseToMat(&sparse, &mat);

    implicitArnoldiSparse(&sparse);

    printf("INPUT:\n");
    printMat(&mat);
    
    /*Sparse sparse;
    Matrix mat, Q, H;
    initSparse(&sparse, 6, 6, 0.2);
    initMat(&mat, sparse.rows, sparse.cols);
    initMat(&Q, sparse.rows, sparse.cols);
    initMat(&H, sparse.rows+1, sparse.cols);
    arnoldiSparse(&sparse, &Q, &H);
    printf("Arnoldi finished.\n");


    sparseToMat(&sparse, &mat);
    printf("Conversion complete.\n");

    printf("sparse:\n");
    printMat(&mat);
    printf("H:\n");
    printMat(&H);
    printf("output:\n");
    qrHess(&H);

    // cleanup
    for (int i = 0; i < H.rows; i++) {
        for (int j = 0; j < H.cols; j++) {
            if (fabs(H.array[i][j]) < 0.000000009)
                H.array[i][j] = 0;
        }
    }

    printMat(&H);

    printf("proposed eigenvalues:\n");
    for (int i = 0, j = 0; i < H.rows-1; i++, j++) {
        double w = H.array[i][i];
        double x = H.array[i][i+1];
        double y = H.array[i+1][i];
        double z = H.array[i+1][i+1];

        double a = 1;
        double b = w + z;
        double c = w*z - y*x;

        double root = b*b - 4*a*c;
        if (root < 0) {
            printf("complex\n");
            printf("eigenvalue %d, part 1: %lf + %lfj\n", j, b/(2*a), sqrtf(-1*root)/(2*a));
            printf("eigenvalue %d, part 2: %lf + %lfj\n", j, b/(2*a), -1*sqrtf(-1*root)/(2*a));
            continue;
        }

        printf("eigenvalue %d, part 1: %lf\n", j, (b + sqrtf(root))/(2*a));
        printf("eigenvalue %d, part 2: %lf\n", j, (b - sqrtf(root))/(2*a));
    }

    if (H.rows % 2 != 0) {
        int i = H.rows - 2;
        double w = H.array[i][i];
        double x = H.array[i][i+1];
        double y = H.array[i+1][i];
        double z = H.array[i+1][i+1];

        double a = 1;
        double b = w + z;
        double c = w*z - y*x;

        double root = b*b - 4*a*c;
        if (root < 0) {
            printf("complex\n");
        }
        else {
            printf("eigenvalue %d, part 1: %lf\n", i, (-1*b + sqrtf(root))/(2*a));
            printf("eigenvalue %d, part 2: %lf\n", i, (-1*b - sqrtf(root))/(2*a));
        }
    }*/
    
    //printf("Spectral Radius: %lf\n", eig);

    /*Matrix mat, vec;
    initMat(&mat, 5, 5);
    initMat(&vec, 5, 1);
    vec.array[0][0] = 1;
    vec.array[1][0] = 1;
    vec.array[2][0] = 1;
    vec.array[3][0] = 1;
    vec.array[4][0] = 1;

    sparseToMat(&sparse, &mat);

    Matrix eye;
    initIdent(&eye, 5, 5);

    set(&eye, 2, 3, 4);

    Matrix out1, out2, out3;
    initMat(&out1, 5, 5);
    initMat(&out2, 5, 5);
    initMat(&out3, 5, 1);

    matDot(&mat, &eye, &out1);
    printf("Check 1. Original vs product:\n");
    printMat(&mat);
    printf("vs\n");
    printMat(&out1);

    sparseDotFirst(&sparse, &eye, &out2);
    printf("Check 2. Original vs product:\n");
    sparsePrint(&sparse);
    printf("vs\n");
    printMat(&out2);

    sparseDotFirst(&sparse, &vec, &out3);
    printf("Check 3. Original vs product:\n");
    printMat(&mat);
    printf("vs\n");
    printMat(&out3);*/
}
