#ifndef ESNFILE
#define ESNFILE
#include <stdlib.h>
#include <stdio.h>
#include <sparse.h>
#include <matrix.h>
#include <linalg.h>
#include <generation.h>

#define SPECMIN 0.001
#define GNUPLOT "gnuplot -persist"

//--------------------------------------------\\
// Net details:                               \\
// - desired sinewave d(n) = 1/2(sin(n/4))    \\
// - no input                                 \\
// - teacher signal 300-step sequence of d(n) \\
// - sigmoid units, activation f = tanh       \\
// - 40 reservoir units                       \\
// - 1 output unit                            \\
//--------------------------------------------\\

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

enum stateUpdateVersion { DEFAULT, STEIL };

int dampening(ESN* esn);
void cleanWeights(Weights* weights, int inputs);

void initWeights(Weights* weights, int inputs, int resSize, int outputs, int batchSize, double alpha) {
    printf("Diagnostic:\n");
    printf("-- inputs: %d\n", inputs);
    printf("-- resSize: %d\n", resSize);
    printf("-- outputs: %d\n", outputs);
    printf("-- batchSize: %d\n", batchSize);
    printf("-- alpha: %lf\n", alpha);

    double density = 1. / (double)resSize;

    weights->reservoir = malloc(sizeof(*(weights->reservoir)));
    weights->outputs = malloc(sizeof(*(weights->outputs)));
    weights->feedback = malloc(sizeof(*(weights->feedback)));

    if (inputs) {
        printf("Inputs initialized.\n");
        weights->inputs = malloc(sizeof(*(weights->inputs)));
        initRandom(weights->inputs, resSize, inputs);
    }
    else
        printf("No inputs.\n");

    int limit = 100;

    initSparse(weights->reservoir, resSize, resSize, density);
    printf("Reservoir initialized.\nScaling...\n");
    // Scale the reservoir
    double spec = spectralRadius(weights->reservoir);
    for (int i = 0; i < limit && fabs(spec) < SPECMIN; i++) {
        printf("Warning: scaling failed (spec == %.15lf). Reinitializing...\n", spec);
        cleanSparse(weights->reservoir);
        initSparse(weights->reservoir, resSize, resSize, density);
        spec = spectralRadius(weights->reservoir);
    }

    //printf("pre-scaled reservoir:\n");
    //sparsePrint(weights->reservoir);
    printf("|spectral radius|: %lf\n", spec);

    scalarSparseDiv(weights->reservoir, spec);
    scalarSparseMult(weights->reservoir, alpha);
    printf("Reservoir scaling complete.\n");

    initIdent(weights->outputs, outputs, resSize + inputs);
    printf("Outputs initialized.\n");

    initRandomNormal(weights->feedback, resSize, outputs);
    //scalarMatMult(weights->feedback, 2.);
    printf("Feedback initialized.\n");
}

void initCollections(Collections* collec, int inputs, int resSize, int outputs, int batchSize) {
    collec->extState = malloc(sizeof(*(collec->extState)));
    collec->extTeacher = malloc(sizeof(*(collec->extTeacher)));

    // these are initialized to one because
    // collected states are appended to the matrix

    // TODO: Consider no appending, and just assignment
    // i.e. initialize these matrices to a fixed size
    // and assign values from there
    initMat(collec->extState, 1, resSize + inputs);
    initMat(collec->extTeacher, 1, outputs);
}

void initStates(States* states, int inputs, int resSize, int outputs) {
    states->currentState = malloc(sizeof(*(states->currentState)));
    states->currentExtState = malloc(sizeof(*(states->currentState)));
    states->currentTeacher = malloc(sizeof(*(states->currentTeacher)));
    states->output = malloc(sizeof(*(states->output)));

    initRandom(states->currentState, 1, resSize);
    initMat(states->currentExtState, 1, resSize + inputs);
    initMat(states->currentTeacher, 1, outputs);
    initMat(states->output, 1, outputs);

    for (int i = 0; i < resSize; i++)
        states->currentExtState->array[0][i] = states->currentState->array[0][i];
}

void initNet(ESN* esn, int inputs, int resSize, int outputs, int batchSize, double alpha, int washout) {
    esn->weights = malloc(sizeof(*(esn->weights)));
    esn->collec = malloc(sizeof(*(esn->collec)));
    esn->states = malloc(sizeof(*(esn->states)));

    initWeights(esn->weights, inputs, resSize, outputs, batchSize, alpha);
    initCollections(esn->collec, inputs, resSize, outputs, batchSize);
    initStates(esn->states, inputs, resSize, outputs);

    esn->inputs = inputs;
    esn->resSize = resSize;
    esn->outputs = outputs;
    esn->batchSize = batchSize;
    esn->washout = washout;

    esn->initialized = 1;

    printf("Looking for echo state property...\n");
    int limit = 500;
    int i = 0;
    for (; i < limit && dampening(esn);) {
        printf("-- Attempt %d\n", ++i);
        cleanWeights(esn->weights, esn->inputs);
        initWeights(esn->weights, inputs, resSize, outputs, batchSize, alpha);
        printf("-- NO ECHO STATE FOUND. Trying again...\n");
    }
    
    printf("Echo state property found.\n");
    printf("Final reservoir:\n");
    sparsePrint(esn->weights->reservoir);

    printf("Initialization complete.\n");
}

void cleanWeights(Weights* weights, int inputs) {
    if (inputs)
        cleanMat(weights->inputs);

    cleanSparse(weights->reservoir);
    cleanMat(weights->outputs);
    cleanMat(weights->feedback);
}

void cleanCollections(Collections* collec) {
    cleanMat(collec->extState);
    cleanMat(collec->extTeacher);
}

void cleanStates(States* states) {
    cleanMat(states->currentState);
    cleanMat(states->currentExtState);
    cleanMat(states->currentTeacher);
    cleanMat(states->output);
}

void cleanNet(ESN* esn) {
    cleanWeights(esn->weights, esn->inputs);
    printf("weights cleaned\n");
    cleanCollections(esn->collec);
    printf("collections cleaned\n");
    cleanStates(esn->states);
    printf("states cleaned\n");
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
    //return sin(3.14*t/8);
    return 0.5*sin(t/4.);
}

// adds a random noise term to each item in the matrix
// -0.001 < x < 0.001
void noise(Matrix* mat) {
    double divisor = 100.;
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++)
            mat->array[i][j] += randomDouble() / divisor;
    }
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

void updateStateNoInput(ESN* esn, int which, int teacherForcing) {
    Matrix res, back;
    States* states = esn->states;
    Weights* weights = esn->weights;

    initMat(&res, states->currentState->rows, weights->reservoir->rows);
    initMat(&back, weights->feedback->rows, states->currentTeacher->size);

    sparseDotFirst(weights->reservoir, states->currentState, &res);
    //printf("res:\n");
    //printMat(&res);
    // change later
    if (teacherForcing) {
        matDot(weights->feedback, states->currentTeacher, &back);
        noise(&back);
    }
    else
        matDot(weights->feedback, states->output, &back);

    if (which == STEIL) {
        double leakRate = 0.05;
        scalarMatMult(states->currentState, 1 - leakRate);
        
        Matrix temp;
        initMat(&temp, res.rows, res.cols);
        matAdd(&res, &back, &temp);
        sigmoid(&temp);

        scalarMatMult(&temp, leakRate);

        matAdd(&temp, states->currentState, states->currentState);

        cleanMat(&temp);
    }

    else {
        zeroMat(states->currentState);

        matAdd(&res, &back, states->currentState);

        sigmoid(states->currentState);
    }

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
        updateStateNoInput(esn, STEIL, 1);
    }

    for (int t = esn->washout; t < esn->batchSize; t++) {
        appendMatRow(collec->extState, states->currentState);
        updateStateNoInput(esn, STEIL, 1);
        sigmoid(states->currentTeacher);
        appendMatRow(collec->extTeacher, states->currentTeacher);

        // TODO: Update for a more general approach (i.e. don't rely on a function)
        for (int i = 0; i < states->currentTeacher->rows; i++) {
            for (int j = 0; j < states->currentTeacher->cols; j++)
                states->currentTeacher->array[i][j] = desired(t);
        }
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

    double alpha = 0.5;

    matDot(&resT, collec->extState, &R);
    //scalarMatDiv(&R, esn->batchSize - esn->washout);
    matDot(&resT, collec->extTeacher, &P);
    //scalarMatDiv(&P, esn->batchSize - esn->washout);

    for (int i = 0; i < R.rows; i++)
        R.array[i][i] += alpha*alpha;

    inverse(&R, &inv);

    matDot(&inv, &P, esn->weights->outputs);

    cleanMat(&resT);
    cleanMat(&R);
    cleanMat(&inv);
    cleanMat(&P);
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
    //printMat(esn->weights->outputs);

    cleanCollections(esn->collec);
    initCollections(esn->collec, esn->inputs, esn->resSize, esn->outputs, esn->batchSize);
}

// no input
double test(ESN* esn) {
    FILE* net;
    FILE* teacher;
    net = fopen("testData.txt", "w");
    teacher = fopen("teacherData.txt", "w");

    if (net == NULL || teacher == NULL) {
        printf("error opening file\n");
        exit(1);
    }

    double sum = 0;
    int testSize = 250;
    for (int i = 0; i < testSize; i++) {
        updateOutput(esn);
        sum += pow((desired(i) - esn->states->output->array[0][0]), 2);
        fprintf(net, "%d %.15lf %.15lf\n", i, desired(i), esn->states->output->array[0][0]);
        fprintf(teacher, "%d %.15lf\n", i, desired(i));
        updateStateNoInput(esn, STEIL, 0);
    }

    fclose(net);
    fclose(teacher);

    double mse = sum / testSize;

    return mse;
}

void sample(ESN* net, int trials, int inputs, int resSize, int outputs, int batchSize, double alpha, int washout) {
    double sigma = 0;
    for (int i = 0; i < trials; i++) {
        initNet(net, inputs, resSize, outputs, batchSize, alpha, washout);
        train(net);
        double mse = test(net);
        sigma += mse;
        printf("Trial %d MSE: %.15lf\n", i+1, mse);
        cleanNet(net);
    }

    double average = sigma / (double)trials;
    printf("Average MSE: %.15lf\n", average);
}

void displayTestData() {
    FILE* gp;
    gp = popen(GNUPLOT, "w");
    fprintf(gp, "plot 'testData.txt' using 1:2 lt rgb \"red\" with lines, 'testData.txt' using 1:3 lt rgb \"blue\" with lines\n");
}

int dampening(ESN* esn) {
    FILE* net;
    net = fopen("dampData.txt", "w");

    if (net == NULL) {
        printf("error opening file\n");
        exit(1);
    }

    Matrix* currentState = esn->states->currentState;
    cleanMat(esn->states->currentState);
    initRandom(currentState, currentState->rows, currentState->cols);
    for (int i = 0; i < 150; i++) {
        updateOutput(esn);
        fprintf(net, "%.15lf\n", esn->states->output->array[0][0]);
        updateStateNoInput(esn, STEIL, 0);
    }

    fclose(net);

    double previous = esn->states->output->array[0][0];
    double convergenceTol = 0.5;
    for (int i = 0; i < 10; i++) {
        updateOutput(esn);
        if (fabs(esn->states->output->array[0][0] - previous) < convergenceTol
         && fabs(esn->states->output->array[0][0]) < convergenceTol)
            return 0; // waveform converged to zero

        updateStateNoInput(esn, STEIL, 0);
    }

    return 1;
}

void displayDampData() {
    FILE* gp;
    gp = popen(GNUPLOT, "w");
    fprintf(gp, "plot 'dampData.txt' using 1 lt rgb \"blue\" with lines\n");
}

#endif
