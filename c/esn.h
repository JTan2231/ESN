#include <stdlib.h>
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

    Matrix* currentState;
    Matrix* currentExtState;
    Matrix* output;
    Collections* collec;
    Weights* weights;
} ESN;

void initWeights(Weights* weights, int inputs, int resSize, int outputs) {
    weights->reservoir = malloc(sizeof weights->reservoir);
    //weights->outputs = malloc(sizeof weights->outputs);
    weights->feedback = malloc(sizeof weights->feedback);

    if (inputs) {
        weights->inputs = malloc(sizeof weights->inputs);
        initRandomNormal(weights->inputs, resSize, inputs);
    }

    initSparse(weights->reservoir, resSize, resSize, 0.15);
    //initRandomNormal(weights->outputs, outputs, resSize + inputs);
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
    if (weights->inputs != NULL) {
        printf("Inputs:\n");
        printMat(weights->inputs);
    }
    else printf("No Inputs.\n");
    printf("Reservoir:\n");
    sparsePrint(weights->reservoir);
    //printf("Outputs:\n");
    //printMat(weights->outputs);
    printf("Feedback:\n");
    printMat(weights->feedback);
}

// applies logistic function to a matrix
void logistic(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            mat->array[i][j] = 1 / (1 + exp(-2 * (mat->array[i][j] - 0.5)));
}

// get the desired output
// for a sin wave generator
double desired(int t) {
    return (1/2)*sin(t/4.);
}

void updateState(ESN* esn, Matrix* nextIn, Matrix* teachOut) {
    Matrix in, res, back;
    initMat(&in, esn->weights->inputs->rows, nextIn->cols);
    initMat(&res, esn->weights->reservoir->rows, esn->currentState->size);
    initMat(&back, esn->weights->feedback->rows, teachOut->cols);

    matDot(esn->weights->inputs, nextIn, &in);
    sparseDotFirst(esn->weights->reservoir, esn->currentState, &res);
    // change later
    matDot(esn->weights->feedback, teachOut, &back);

    zeroMat(esn->currentState);

    matAdd(&in, &res, esn->currentState);
    matAdd(esn->currentState, &back, esn->currentState);
    
    logistic(esn->currentState);
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
