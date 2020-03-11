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
    int step;
    Vector currentState;
    Vector currentExtState;
    Vector output;
    Matrix extStateCollection;
    Matrix extTeacherCollection;
    Weights* weights;
} ESN;

void initWeights(Weights* weights, int inputs, int resSize, int outputs) {
    weights->inputs = malloc(sizeof weights->inputs);
    weights->reservoir = malloc(sizeof weights->reservoir);
    weights->outputs = malloc(sizeof weights->outputs);
    weights->feedback = malloc(sizeof weights->feedback);
}

// TODO: 0 < spectral radius < 1
void initNet(ESN* esn, int inputs, int resSize, int outputs) {
    esn->weights = malloc(sizeof esn->weights);
    initWeights(esn->weights, inputs, resSize, outputs);

    initRandomNormal(esn->weights->inputs, resSize, inputs);
    initSparse(esn->weights->reservoir, resSize, resSize, 0.1);
    initRandomNormal(esn->weights->outputs, outputs, resSize + inputs);
    initRandomNormal(esn->weights->feedback, resSize, outputs);
    printf("Initialized\n");
}

void printWeights(ESN* esn) {
    printf("Inputs:\n");
    printMat(esn->weights->inputs);
    printf("Reservoir:\n");
    sparsePrint(esn->weights->reservoir);
    printf("Outputs:\n");
    printMat(esn->weights->outputs);
    printf("Feedback:\n");
    printMat(esn->weights->feedback);
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
