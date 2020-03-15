#include <stdio.h>
#include "map.h"
#include "matrix.h"

typedef struct {
    Matrix* inputs;
    Sparse* reservoir;
    Matrix* outputs;
    Matrix* feedback;
} Weights;

typedef struct {
    Matrix* extState;
    Matrix* extTeacher;
} Collections;

typedef struct {
    int step;
    int batchSize;
    int inputs;
    int resSize;
    int outputs;

    Vector currentState;
    Vector currentExtState;
    Vector output;
    Collections* collec;
    Weights* weights;
} ESN;

void initWeights(Weights* weights, int inputs, int resSize, int outputs) {
    weights->inputs = malloc(sizeof weights->inputs);
    weights->reservoir = malloc(sizeof weights->reservoir);
    weights->outputs = malloc(sizeof weights->outputs);
    weights->feedback = malloc(sizeof weights->feedback);

    initRandomNormal(weights->inputs, resSize, inputs);
    initSparse(weights->reservoir, resSize, resSize, 0.35);
    initRandomNormal(weights->outputs, outputs, resSize + inputs);
    initRandomNormal(weights->feedback, resSize, outputs);
}

void initCollections(Collections* collec, int inputs, int resSize, int outputs, int batchSize) {
    collec->extState = malloc(sizeof collec->extState);
    collec->extTeacher = malloc(sizeof collec->extTeacher);

    initMat(collec->extState, batchSize, resSize + inputs);
    initMat(collec->extTeacher, batchSize, outputs);
}

// TODO: 0 < spectral radius < 1
void initNet(ESN* esn, int inputs, int resSize, int outputs, int batchSize) {
    esn->weights = malloc(sizeof esn->weights);
    esn->collec = malloc(sizeof esn->collec);

    initWeights(esn->weights, inputs, resSize, outputs);
    initCollections(esn->collec, inputs, resSize, outputs, batchSize);

    esn->inputs = inputs;
    esn->resSize = resSize;
    esn->outputs = outputs;
    esn->batchSize = batchSize;

    printf("Initialized\n");
}

void printWeights(Weights* weights) {
    printf("Inputs:\n");
    printMat(weights->inputs);
    printf("Reservoir:\n");
    sparsePrint(weights->reservoir);
    printf("Outputs:\n");
    printMat(weights->outputs);
    printf("Feedback:\n");
    printMat(weights->feedback);
}

// applies logistic function to a vector
void logistic(Vector* vec) {
    for (int i = 0; i < vec->size; i++)
        vec->array[i] = 1 / (1 + exp(-2 * (vec->array[i] - 0.5)));
}

// updates output weights using ridge regression
//void updateWeights(ESN* esn)

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
