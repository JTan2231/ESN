#ifndef ESNFILE
#define ESNFILE
#include <stdlib.h>
#include <stdio.h>
#include "map.h"
#include "matrix.h"
#include "linalg.h"
#include "generation.h"

#define SPECTOL 0.001

//--------------------------------------------\\
// Net details:                               \\
// - desired sinewave d(n) = 1/2(sin(n/4))    \\
// - no input                                 \\
// - teacher signal 300-step sequence of d(n) \\
// - sigmoid units, activation f = tanh       \\
// - 40 reservoir units                       \\
// - 1 output unit                            \\
//--------------------------------------------\\

// TODO: Extended system states; accommodate inputs

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
    Matrix* currentState;
    Matrix* currentExtState;
    Matrix* currentTeacher;
    Matrix* output;
} States;

typedef struct {
    int initialized;
    int step;
    int batchSize;
    int washout;
    int inputs;
    int resSize;
    int outputs;

    States* states;
    Collections* collec;
    Weights* weights;
} ESN;

void initWeights(Weights* weights, int inputs, int resSize, int outputs, int batchSize, double alpha) {
    printf("Diagnostic:\n");
    printf("-- inputs: %d\n", inputs);
    printf("-- resSize: %d\n", resSize);
    printf("-- outputs: %d\n", outputs);
    printf("-- batchSize: %d\n", batchSize);
    printf("-- alpha: %lf\n", alpha);

    double density = 0.05;

    weights->reservoir = malloc(sizeof(weights->reservoir));
    weights->outputs = malloc(sizeof(weights->outputs));
    weights->feedback = malloc(sizeof(weights->feedback));

    if (inputs) {
        printf("Inputs initialized.\n");
        weights->inputs = malloc(sizeof(weights->inputs));
        initRandom(weights->inputs, resSize, inputs);
    }
    else
        printf("No inputs.\n");

    int limit = 100;

    initSparse(weights->reservoir, resSize, resSize, density);
    printf("Reservoir initialized.\nScaling...\n");
    // Scale the reservoir
    double spec = fabs(spectralRadius(weights->reservoir));
    for (int i = 0; i < limit && spec < SPECTOL; i++) {
        printf("Warning: scaling failed. Reinitializing...\n");
        cleanSparse(weights->reservoir);
        initSparse(weights->reservoir, resSize, resSize, density);
        spec = spectralRadius(weights->reservoir);
    }

    printf("pre-scaled reservoir:\n");
    sparsePrint(weights->reservoir);
    printf("spectral radius: %lf\n", spec);

    scalarSparseDiv(weights->reservoir, spec);
    scalarSparseMult(weights->reservoir, alpha);
    printf("Reservoir scaling complete.\n");

    initRandom(weights->outputs, outputs, resSize + inputs);
    printf("Outputs initialized.\n");

    initRandom(weights->feedback, resSize, outputs);
    printf("Feedback initialized.\n");
}

void initCollections(Collections* collec, int inputs, int resSize, int outputs, int batchSize) {
    collec->extState = malloc(sizeof(collec->extState));
    collec->extTeacher = malloc(sizeof(collec->extTeacher));

    // these are initialized to one because
    // collected states are appended to the matrix

    // TODO: Consider no appending, and just assignment
    // i.e. initialize these matrices to a fixed size
    // and assign values from there
    initMat(collec->extState, 1, resSize + inputs);
    initMat(collec->extTeacher, 1, outputs);
}

void initStates(States* states, int inputs, int resSize, int outputs) {
    states->currentState = malloc(sizeof(states->currentState));
    states->currentExtState = malloc(sizeof(states->currentState));
    states->currentTeacher = malloc(sizeof(states->currentTeacher));
    states->output = malloc(sizeof(states->output));

    initRandom(states->currentState, 1, resSize);
    initMat(states->currentExtState, 1, resSize + inputs);
    initMat(states->currentTeacher, 1, outputs);
    initMat(states->output, 1, outputs);

    for (int i = 0; i < resSize; i++)
        states->currentExtState->array[0][i] = states->currentState->array[0][i];
}

void initNet(ESN* esn, int inputs, int resSize, int outputs, int batchSize, double alpha, int washout) {
    esn->weights = malloc(sizeof(esn->weights));
    esn->collec = malloc(sizeof(esn->collec));
    esn->states = malloc(sizeof(esn->states));

    initWeights(esn->weights, inputs, resSize, outputs, batchSize, alpha);
    initCollections(esn->collec, inputs, resSize, outputs, batchSize);
    initStates(esn->states, inputs, resSize, outputs);

    esn->inputs = inputs;
    esn->resSize = resSize;
    esn->outputs = outputs;
    esn->batchSize = batchSize;
    esn->washout = washout;

    esn->initialized = 1;

    printf("Initialization complete.\n");
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
}

void sigmoid(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            mat->array[i][j] = tanh(mat->array[i][j]);
    }
}

// applies arctanh(x) to a matrix
void inverseSigmoid(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            mat->array[i][j] = atanh(mat->array[i][j]);
    }
}

// get the desired output
// for a sin wave generator
double desired(int t) {
    return (1./2.)*sin(t/4.);
}

void updateState(ESN* esn, Matrix* nextIn) {
    Matrix in, res, back;
    
    States* states = esn->states;
    Weights* weights = esn->weights;

    initMat(&in, weights->inputs->rows, nextIn->cols);
    initMat(&res, weights->reservoir->rows, states->currentState->rows);
    initMat(&back, weights->feedback->rows, states->currentTeacher->size);

    matDot(weights->inputs, nextIn, &in);
    sparseDotFirst(weights->reservoir, states->currentState, &res);
    // change later
    matDot(weights->feedback, states->currentTeacher, &back);

    zeroMat(states->currentState);

    matAdd(&in, &res, states->currentState);
    matAdd(states->currentState, &back, states->currentState);

    logistic(states->currentState);

    cleanMat(&in);
    cleanMat(&res);
    cleanMat(&back);
}

void updateStateNoInput(ESN* esn, int teacherForcing) {
    Matrix res, back;

    States* states = esn->states;
    Weights* weights = esn->weights;

    initMat(&res, states->currentState->rows, weights->reservoir->rows);
    initMat(&back, weights->feedback->rows, states->currentTeacher->size);

    sparseDotFirst(weights->reservoir, states->currentState, &res);
    //printf("res:\n");
    //printMat(&res);
    // change later
    if (teacherForcing)
        matDot(weights->feedback, states->currentTeacher, &back);
    else
        matDot(weights->feedback, states->output, &back);

    //printf("back:\n");
    //printMat(&back);

    zeroMat(states->currentState);

    matAdd(&res, &back, states->currentState);

    sigmoid(states->currentState);

    cleanMat(&res);
    cleanMat(&back);
}

// no input
void updateOutput(ESN* esn) {
    Weights* weights = esn->weights;
    States* states = esn->states;

    matDot(weights->outputs, states->currentState, states->output);
}

void collectStates(ESN* esn) {
    States* states = esn->states;
    Weights* weights = esn->weights;
    Collections* collec = esn->collec;

    for (int t = 0; t < esn->washout; t++) {
        states->currentTeacher->array[0][0] = desired(t);
        updateStateNoInput(esn, 1);
    }

    for (int t = esn->washout; t < esn->batchSize; t++) {
        /*printf("currentTeacher:\n");
        printMat(states->currentTeacher);
        printf("currentState:\n");
        printMat(states->currentState);*/

        appendMatRow(collec->extState, states->currentState);
        updateStateNoInput(esn, 1);
        inverseSigmoid(states->currentTeacher);
        appendMatRow(collec->extTeacher, states->currentTeacher);

        Matrix newTeach;
        initMat(&newTeach, states->currentTeacher->rows, states->currentTeacher->cols);
        // TODO: Update for a more general approach (i.e. don't rely on a function)
        for (int i = 0; i < newTeach.rows; i++) {
            for (int j = 0; j < newTeach.cols; j++)
                newTeach.array[i][j] = desired(t);
        }

        cleanMat(states->currentTeacher);
        states->currentTeacher->array = newTeach.array;
    }
}

// uses Tikhonov Regularization to calculate output weights
void tikhonov(ESN* esn) {
    Collections* collec = esn->collec;

    Matrix resT, R, inv, P;
    initTranspose(collec->extState, &resT);
    initMat(&R, resT.rows, collec->extState->cols);
    initMat(&inv, R.rows, R.cols);
    initMat(&P, resT.rows, collec->extTeacher->cols);

    double alpha = 1.3;

    matDot(&resT, collec->extState, &R);
    //scalarMatDiv(&R, esn->batchSize - esn->washout);
    matDot(&resT, collec->extTeacher, &P);
    //scalarMatDiv(&P, esn->batchSize - esn->washout);

    for (int i = 0; i < R.rows; i++)
        R.array[i][i] += alpha*alpha;

    inverse(&R, &inv);

    matDot(&inv, &P, esn->weights->outputs);
}

void train(ESN* esn) {
    assert(esn->initialized);
    printf("Training:\n");

    printf("-- Collecting states...\n");
    collectStates(esn);
    printf("-- Success: States collected.\n");
    //printMat(esn->collec->extState);

    printf("-- Calculating output weights...\n");
    tikhonov(esn);
    printf("-- Success: Output weights calculated.\n");
    printMat(esn->weights->outputs);
}

// no input
void test(ESN* esn) {
    for (int i = 0; i < esn->states->currentTeacher->rows; i++) {
        for (int j = 0; j < esn->states->currentTeacher->cols; j++)
            esn->states->output->array[i][j] = esn->states->currentTeacher->array[i][j];
    }

    double sum = 0;
    for (int i = 0; i < 500; i++) {
        updateOutput(esn);
        sum += pow((desired(i+esn->batchSize) - esn->states->output->array[0][0]), 2);
        printf("d(n), y(n): %lf, %lf\n", desired(i+300), esn->states->output->array[0][0]);
        //printf("Current state:\n");
        //printMat(esn->states->currentState);
        updateStateNoInput(esn, 0);
    }

    double mse = sum / 50.;

    printf("Error: %.15lf\n", mse);
}

void dampening(ESN* esn) {
    Matrix* currentState = esn->states->currentState;
    cleanMat(esn->states->currentState);
    initRandom(currentState, currentState->rows, currentState->cols);
    for (int i = 0; i < 50; i++) {
        updateOutput(esn);
        printf("output %d: %.15lf\n", i, esn->states->output->array[0][0]);
        updateStateNoInput(esn, 0);
    }
}

#endif
