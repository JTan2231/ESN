#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "matrix.h"

int main() {
    srand(time(0));

    Matrix m1;
    initRandomNormal(&m1, 7, 7);

    printMat(&m1);
    printf("\n");
}
