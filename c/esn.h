#include <stdio.h>
#include "map.h"
#include "matrix.h"

typedef struct {
    int step;
    Vector currentState;
    Vector currentExtState;
    Vector output;
    Matrix extStateCollection;
    Matrix extTeacherCollection;
    Matrix inputs;
    Sparse reservoir;
    Matrix outputs;
    Matrix feedback;
} ESN;

// TODO: 0 < spectral radius < 1
void initWeights(ESN* esn, int inputs, int resSize, int outputs) {
    int rows = 10;
    int cols = 10;
    printf("Initialized\n");
    initRandomNormal(&(esn->inputs), resSize, inputs);
    initSparse(&(esn->reservoir), resSize, resSize, 0.1);
    initRandomNormal(&(esn->outputs), outputs, resSize + inputs);
    initRandomNormal(&(esn->feedback), resSize, outputs);
}

void printWeights(ESN* esn) {
    printf("Inputs:\n");
    printMat(&(esn->inputs));
    printf("Reservoir:\n");
    sparsePrint(&(esn->reservoir));
    printf("Outputs:\n");
    printMat(&(esn->outputs));
    printf("Feedback:\n");
    printMat(&(esn->feedback));
}

// applies logistic function to a vector
void logistic(Vector* vec) {
    for (int i = 0; i < vec->size; i++)
        vec->array[i] = 1 / (1 + exp(-2 * (vec->array[i] - 0.5)));
}

/*void updateState(ESN* esn, Vector* nextInput) {
    Vector a, b, c, d, e;
    initVec(&a, esn->currentState->size);
    initVec(&b, nextInput->size);
    initVec(&c, esn->output->size);
    initVec(&d, a.size);

    sparseVecDot(esn->weights->reservoir, esn->currentState, &a);
    vecDot(esn->weights->inputs,
*/
/*void computeOutput(ESN* esn) {
    matVecDot(esn->outputs, esn->currentExtState, esn->output);
}*/
