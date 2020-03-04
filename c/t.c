#include <stdio.h>
#include <time.h>
#include "matrix.h"
#include "generation.h"

int main() {
    srand(time(0));

    Matrix mat;
    Matrix inv;
    Matrix lower;
    Matrix upper;
    Matrix check;

    initRandom(&mat, 5, 5);
    initMat(&inv, 5, 5);
    initMat(&lower, 5, 5);
    initMat(&upper, 5, 5);
    initMat(&check, 5, 5);

    crout(&mat, &lower, &upper);
    printf("Check\n");
    inverse(&mat, &lower, &upper, &inv);
    matDot(&mat, &inv, &check);

    /*printf("Input:\n");
    printMat(&mat);
    printf("Upper:\n");
    printMat(&upper);
    printf("Lower:\n");
    printMat(&lower)*/;
    printf("Inverse:\n");
    printMat(&inv);
    printf("Check:\n");
    printMat(&check);
}
