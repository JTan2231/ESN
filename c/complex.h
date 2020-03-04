#ifndef COMPLEX
#define COMPLEX 
#include <stdio.h>

typedef struct {
    double real;
    double i;
} Complex;

Complex initComp(double real, double i) {
    Complex output;
    output.real = real;
    output.i = i;

    return output;
}

Complex conjugate(Complex comp) {
    Complex output;
    output.real = comp.real;
    output.i = -1*comp.i;

    return output;
}

Complex addScalarComp(double scalar, Complex comp) {
    Complex output;
    output.real = comp.real + scalar;
    output.i = comp.i;

    return output;
}

Complex addComp(Complex c1, Complex c2) {
    Complex output;
    output.real = c1.real + c2.real;
    output.i = c1.i + c2.i;

    return output;
}

Complex multScalarComp(double scalar, Complex comp) {
    Complex output;
    output.real = scalar * comp.real;
    output.i = scalar * comp.i;

    return output;
}

Complex multComp(Complex c1, Complex c2) {
    Complex output;
    output.real = c1.real * c2.real - c1.i * c2.i;
    output.i = c1.real * c2.real + c1.i * c2.i;

    return output;
}

Complex divComp(Complex c1, Complex c2) {
    Complex output;
    double den = c1.real * c1.real + c2.real * c2.real;
    assert(den != 0);
    
    output.real = (c1.real * c2.real + c1.i * c2.i) / den;
    output.i = (c1.real * c2.i - c2.real * c1.i) / den;

    return output;
}

Complex recipComp(Complex comp) {
    Complex output;
    double den = comp.real * comp.real + comp.i * comp.i;
    assert(den != 0);

    output.real = comp.real / den;
    output.i = -1*(comp.i / den);

    return output;
}

#endif
