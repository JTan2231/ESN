#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "matrix.h"
#include "map.h"
#include "generation.h"

int main() {
    srand(time(0));
    
    /*Parent p;
    initParent(&p, 5);

    for (int i = 0; i < 5; i++) {
        printf("Before: \n");
        parentPrint(&p);
        parentAdd(&p, randRange(10), marsagliaPolar());
        parentAdd(&p, randRange(10), marsagliaPolar());
        printf("After: \n");
        parentPrint(&p);
    }*/

    ParentArray mat;
    initSparse(&mat, 10, 10, 0.1);
    Matrix vec;
    initRandom(&vec, 1, 10);
    Matrix out;
    initMat(&out, 1, 10);

    printf("Sparse: \n");
    parentArrayPrint(&mat);
    printf("vec: \n");
    printMat(&vec);

    /*sparseDotFirst(&mat, &vec, &out);

    printf("out: \n");
    printMat(&out);*/

    eigenVector(&mat, &vec, 10);

    printf("The eigenvector of \n");
    parentArrayPrint(&mat);
    printf("is \n");
    printMat(&vec);
    

    /*ParentArray pa;
    initSparse(&pa, 10, 10, 0.1);
    Matrix m;
    Matrix out;
    initRandomNormal(&m, 10, 10);
    initMat(&out, 10, 10);
    parentArrayPrint(&pa);
    printMat(&m);
    sparseDotSecond(&m, &pa, &out);
    printf("out:\n");
    printMat(&out);*/
    
    /*
    for (int i = 0; i < 10; i++) {
        printf("Before: \n");
        parentArrayPrint(&pa);
        Parent p;
        initParent(&p, 5);
        parentAdd(&p, randRange(10), marsagliaPolar());
        parentAdd(&p, randRange(10), marsagliaPolar());
        parentArrayAdd(&pa, &p);
        printf("After: \n");
        parentArrayPrint(&pa);
    }*/
    
    return 0;
}
