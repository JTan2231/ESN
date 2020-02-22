#include <stdlib.h>
#include <math.h>

const float PI = 3.14159;
const float RAND_UPPER = 3;

// basically a coin toss
int randBinary() {
    const int divisor = RAND_MAX/(2);
    int retval;

    do {
        retval = rand() / divisor;
    } while (retval > 1);

    return retval;
}

float randomFloat() {
    float output = (float)rand()/(float)(RAND_MAX/RAND_UPPER);
    if (randBinary()) return -1*output;
    return output;
}

float normal(float x) {
    const float mean = 0;
    const float variance = 1;

    return expf(-0.5 * pow((x-mean)/variance, 2)) / (variance * sqrtf(2*PI));
}

// generate value from a normal distribution
// see Marsaglia Polar Method
float marsagliaNormal() {
    float U, V, S;
    do {
        U = (float)rand()/(float)(RAND_MAX);
        if (randBinary()) U *= -1;
        V = (float)rand()/(float)(RAND_MAX);
        if (randBinary()) V *= -1;
        S = U*U + V*V;
    } while (S >= 1);

    return U * sqrtf(-2*logf(S)/S);
}

