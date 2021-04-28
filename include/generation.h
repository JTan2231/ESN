#ifndef GENERATION
#define GENERATION

//------------------------------------\\
// Random Number Generation           \\
//------------------------------------\\

// basically a coin toss
int randBinary();

int randRange(int max);

double randomDouble();

double randomPositiveDouble();

//------------------------------------\\
// Statistics                         \\
//------------------------------------\\

// TODO: -1 < x < 1
double normal(double x);

// generate value from a normal distribution
// see Marsaglia Polar Method
double marsagliaPolar();

#endif
