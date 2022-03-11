#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <matrix.h>
#include <sparse.h>
#include <generation.h>
#include <linalg.h>
#include <esn.h>

int main() {
    srand(time(0));

    ESN net;

    int inputs = 0;
    int resSize = 20;
    int outputs = 1;
    int batchSize = 300;
    int washout = 100;
    double alpha = 20.;

    initNet(&net, inputs, resSize, outputs, batchSize, alpha, washout);
    train(&net);
    double mse = test(&net);
    displayTestData();

    printf("mse: %.15lf\n", mse);
    
    return 0;
}
