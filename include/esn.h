#ifndef ESNFILE
#define ESNFILE
#include <stdlib.h>
#include <stdio.h>
#include "sparse.h"
#include "matrix.h"
#include "linalg.h"
#include "generation.h"

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

enum stateUpdateVersion { DEFAULT, STEIL };

int dampening(ESN* esn);
void cleanWeights(Weights* weights, int inputs);

//--------------------------------------------------------\\
// Initialization Functions                               \\
//--------------------------------------------------------\\
// These functions initialize memory buffers and their    \\
// initial values using the given parameters.             \\
//--------------------------------------------------------\\

void initWeights(Weights* weights, int inputs, int resSize, int outputs, int batchSize, double alpha);
void initCollections(Collections* collec, int inputs, int resSize, int outputs, int batchSize);
void initStates(States* states, int inputs, int resSize, int outputs);
void initNet(ESN* esn, int inputs, int resSize, int outputs, int batchSize, double alpha, int washout);

//--------------------------------------------------------\\
// Cleanup Functions                                      \\
//--------------------------------------------------------\\
// Cleanup functions corresponding to the initializers    \\
// above. Deallocates allocated pointers.                 \\
//--------------------------------------------------------\\

void cleanWeights(Weights* weights, int inputs);
void cleanCollections(Collections* collec);
void cleanStates(States* states);
void cleanNet(ESN* esn);

void printWeights(Weights* weights);

//--------------------------------------------------------\\
// Math Functions                                         \\
//--------------------------------------------------------\\
// Applies element-wise mathematical operations to a      \\
// matrix.                                                \\
//--------------------------------------------------------\\

// Logistic function, defined as f(x;L,k,y) = L / -(1 + exp(-k*(x-y))
// where:
// -- L = 1
// -- k = -2
// -- y = 0.5
void logistic(Matrix* mat);

// Sigmoid function, defined as f(x) = 1 / (1 + exp(-x))
void sigmoid(Matrix* mat);

// Inverse hyperbolic tangent, defined as f(x) = arctanh(x)
void inverseSigmoid(Matrix* mat);

// Sine function, defined here as f(t;k,l) = k*sin(PI*t/l)
double desired(int t);

// Adds noise added from a uniform distribution, -0.001 < x < 0.001
void noise(Matrix* mat);

//--------------------------------------------------------\\
// Update Functions                                       \\
//--------------------------------------------------------\\
// Various functions for updating the inner states of     \\
// the echo state network. States in this context means   \\
// the weight matrices
//--------------------------------------------------------\\

// Updates the reservoir state using function (1) defined in 
// http://www.scholarpedia.org/article/Echo_state_network#Formalism_and_theory
// with an external input, Matrix* nextIn
void updateState(ESN* esn, Matrix* nextIn);

// Updates the reservoir state using the same function as above,
// though this time without an external input. Used during training.
//
// Parameters:
// -- int which: method of state update
// -- int teacherForcing: boolean for using teacher forcing
void updateStateNoInput(ESN* esn, int which, int teacherForcing);

// Calculates the network output using the current reservoir state
// and output weights
void updateOutput(ESN* esn);

// Function for conditioning the reservoir state during training
// as described in
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.378.4095&rep=rep1&type=pdf
void collectStates(ESN* esn);

// Tikhonov Regularization to calculate output weights.
// Used at the end of training.
void tikhonov(ESN* esn);

// Wrapper function to ease the training process.
// Trains the given echo state network.
void train(ESN* esn);

// Test a trained (or untrained) echo state network
// on sine function approximation and write the results
// to file.
double test(ESN* esn);

// Train a group of echo state networks and view the results.
void sample(ESN* net, int trials, int inputs, int resSize, int outputs, int batchSize, double alpha, int washout);

// Graphical display of test data.
void displayTestData();

// Check if the echo state network has the damping property.
int dampening(ESN* esn);
void displayDampData();

#endif
