/*
 * Author: Joey Tan
 * Date Created: 2-19-20
 * Last Edit: 2-23-20, Joey Tan
 */

#ifndef GENERATION
#define GENERATION
#include <stdlib.h>
#include <math.h>

const double PI = 3.14159;
const double RAND_UPPER = 1;

//------------------------------------\\
// Random Number Generation           \\
//------------------------------------\\

// basically a coin toss
int randBinary() {
    const int divisor = RAND_MAX/(2);
    int retval;

    do {
        retval = rand() / divisor;
    } while (retval > 1);

    return retval;
}

int randRange(int max) {
    const int divisor = RAND_MAX/(max);
    int retval;

    do {
        retval = rand() / divisor;
    } while (retval > max);

    return retval;
}

double randomDouble() {
    double output = (double)rand()/(double)(RAND_MAX/RAND_UPPER);
    if (randBinary()) return -1*output;
    return output;
}

//------------------------------------\\
// Statistics                         \\
//------------------------------------\\

// TODO: -1 < x < 1
double normal(double x) {
    const double mean = 0;
    const double variance = 1;

    return expf(-0.5 * pow((x-mean)/variance, 2)) / (variance * sqrtf(2*PI));
}

// generate value from a normal distribution
// see Marsaglia Polar Method
double marsagliaPolar() {
    double U, V, S;
    do {
        U = (double)rand()/(double)(RAND_MAX);
        if (randBinary()) U *= -1;
        V = (double)rand()/(double)(RAND_MAX);
        if (randBinary()) V *= -1;
        S = U*U + V*V;
    } while (S >= 1);

    return U * sqrtf(-2*logf(S)/S);
}

#endif
