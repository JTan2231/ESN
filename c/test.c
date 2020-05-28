#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "matrix.h"
#include "sparse.h"
#include "generation.h"
#include "linalg.h"
#include "esn.h"

// rudimentary test file
// to make sure everything works as expected
// unfinished

// TODO: Function testing
// TODO: Expected outputs from given inputs

int main() {
    srand(time(0));

    ESN net;

    int inputs = 0;
    int resSize = 10;
    int outputs = 1;
    int batchSize = 500;
    int washout = 50;
    double alpha = 0.75;

    initNet(&net, inputs, resSize, outputs, batchSize, alpha, washout);
    train(&net);
    double mse = test(&net);
    displayTestData();

    printf("mse: %.15lf\n", mse);
    
    return 0;
}
