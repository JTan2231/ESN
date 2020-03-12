#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "matrix.h"
#include "map.h"
#include "generation.h"
#include "linalg.h"
#include "esn.h"

// rudimentary test file
// to make sure everything works as expected
// unfinished

// TODO: Function testing
// TODO: Expected outputs from given inputs

void printHeader(char header[]) {
    printf("-----------------\n%s\n-----------------\n", header);
}

void printSubHeader(char subHeader[]) {
    printf("\n>>> %s\n\n", subHeader);
}

int main() {
    srand(time(0));

    /*ESN esn;
    initNet(&esn, 1, 25, 1);
    //printWeights(&esn);

    Matrix m1, m2, m3;
    initRandomNormal(&m1, 5, 7);
    initRandomNormal(&m2, 7, 1);
    initMat(&m3, 1, 5);

    matDot(&m1, &m2, &m3);
    printMat(&m3);*/

    Matrix mat;
    initMat(&mat, 5, 5);
    //initMat(&mat, 5, 5);
    //set(&mat, -0.17, 0, 2);
    //set(&mat, 2.26, 4, 2);
    Matrix arnD, arnDH;
    initMat(&arnD, mat.rows, mat.cols);
    initMat(&arnDH, mat.rows+1, mat.cols);

    Sparse sparse;
    initSparse(&sparse, 5, 5, 0.35);
    Matrix arnSQ, arnSH;
    initMat(&arnSQ, sparse.rows, sparse.cols);
    initMat(&arnSH, sparse.rows+1, sparse.cols);
    
    Matrix a, b;
    initMat(&a, 5, 5);
    initMat(&b, 5, 5);
    
    arnoldiSparse(&sparse, &arnSQ, &arnSH);
    
    printf("Sparse:\n");
    sparsePrint(&sparse);
    printf("Hessenberg:\n");
    printMat(&arnSH);
    /*sparseToMat(&sparse, &mat);
    printf("Converted Sparse:\n");
    printMat(&mat);*/
    
    /*shrinkMat(&arnSH, arnSQ.rows, arnSQ.cols);
    
    printf("Sparse x Q:\n");
    matDot(&mat, &arnSQ, &a);
    printMat(&a);
    printf("Q x H:\n");
    matDot(&arnSQ, &arnSH, &b);
    printMat(&b);*/

    /*arnoldiDense(&mat, &arnD, &arnDH);
    printf("Dense Orthogonal:\n");
    printMat(&arnD);
    printf("Dense Hessenberg:\n");
    printMat(&arnDH);*/

    Matrix test;
    /*printf("Dense Orthogonality Check:\n");
    Matrix arnDT;
    initTranspose(&arnD, &arnDT);
    initMat(&test, 50, 50);
    matDot(&arnD, &arnDT, &test);
    printMat(&test);
    cleanMat(&test);*/

    /*arnoldiSparse(&sparse, &arnSQ, &arnSH);
    Matrix arnSQT;
    initTranspose(&arnSQ, &arnSQT);
    initMat(&test, arnSQ.rows, arnSQ.cols);
    printf("Sparse:\n");
    sparsePrint(&sparse);*/
    /*printf("ArnoldiSparse:\n");
    printMat(&arnSQ);
    printf("Sparse Orthogonal:\n");
    printMat(&arnSQ);*/
    
    /*Matrix temp, out;
    initMat(&temp, sparse.rows, sparse.cols);
    sparseDotSecond(&arnSQT, &sparse, &temp);
    initMat(&out, sparse.rows, sparse.cols);
    matDot(&temp, &arnSQ, &out);
    printf("Sparse Hessenberg:\n");
    printMat(&arnSH);*/
    
    /*printf("Transpose:\n");
    printMat(&arnSQT);*/
    
    //shrinkMat(&arnSH, arnSQ.rows, arnSQ.cols);
    
    /*printf("Shrunk Hessenberg:\n");
    printMat(&arnSH);
    
    printf("Test:\n");
    printMat(&out);*/
    
    /*printf("QR Algorithm:\n");
    qrSparse(&arnSH, &arnSQ, &test);
    printMat(&test);*/

    /*printf("Sparse Orthogonality Check 1:\n");
    //printf("Transpose:\n");
    //printMat(&arnSQT);
    initMat(&test, arnSQT.rows, arnSQ.cols);
    matDot(&arnSQT, &arnSQ, &test);
    printMat(&test);
    
    cleanMat(&test);
    
    printf("Sparse Orthogonality Check 2:\n");
    initMat(&test, arnSQ.rows, arnSQT.cols);
    matDot(&arnSQ, &arnSQT, &test);
    printMat(&test);*/

    //printHeader("Matrix Initializations");
    /*Matrix mat;
    Matrix mat2;
    Matrix out;
    Matrix transpose;
    ParentArray sparse;*/

    /*printSubHeader("Zeros Matrix Initialization");
    initMat(&mat, 5, 5);
    printMat(&mat);
    cleanMat(&mat);

    printSubHeader("Random Float Matrix Initialization");
    initRandom(&mat, 5, 5);
    printMat(&mat);
    cleanMat(&mat);

    printSubHeader("Random Range Matrix Initialization");
    initRandomRange(&mat, 5, 5, 50);
    printMat(&mat);
    cleanMat(&mat);

    printSubHeader("Random Normal Matrix Initialization");
    initRandomNormal(&mat, 5, 5);
    printMat(&mat);
    cleanMat(&mat);

    printSubHeader("Sparse Matrix Initialization");
    initSparse(&sparse, 7, 7, 0.11);
    parentArrayPrint(&sparse);
    cleanParentArray(&sparse);

    printSubHeader("Matrix Transpose Initialization");
    initRandomNormal(&mat, 5, 7);
    initTranspose(&mat, &transpose);
    
    printf("Matrix:\n");
    printMat(&mat);
    printf("Tranpose:\n");
    printMat(&transpose);

    cleanMat(&mat);
    cleanMat(&transpose);*/

    /*printHeader("Matrix Operations");

    printSubHeader("Square Matrix Multiplication");
    initRandomRange(&mat, 3, 3, 5);
    initRandomRange(&mat2, 3, 3, 5);
    initMat(&out, 3, 3);
    matDot(&mat, &mat2, &out);

    printf("Matrix 1:\n");
    printMat(&mat);
    printf("Matrix 2:\n");
    printMat(&mat2);
    printf("Output:\n");
    printMat(&out);

    cleanMat(&mat);
    cleanMat(&mat2);
    cleanMat(&out);

    printSubHeader("Rectangular Matrix Multiplication");
    initRandomRange(&mat, 3, 5, 5);
    initRandomRange(&mat2, 5, 7, 5);
    initMat(&out, 3, 7);
    matDot(&mat, &mat2, &out);

    printf("Matrix 1:\n");
    printMat(&mat);
    printf("Matrix 2:\n");
    printMat(&mat2);
    printf("Output:\n");
    printMat(&out);

    cleanMat(&mat);
    cleanMat(&mat2);
    cleanMat(&out);

    printSubHeader("Sparse Matrix by Dense Matrix Multiplication");
    initSparse(&sparse, 5, 6, 0.25);
    initRandomNormal(&mat, 3, 5);
    initMat(&out, 3, 6);
    sparseDotSecond(&mat, &sparse, &out);

    printf("Dense Matrix:\n");
    printMat(&mat);
    printf("Sparse Matrix:\n");
    parentArrayPrint(&sparse);
    printf("Output:\n");
    printMat(&out);

    cleanMat(&mat);
    cleanParentArray(&sparse);
    cleanMat(&out);*/

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

    /*ParentArray mat;
    initSparse(&mat, 25, 25, 0.1);
    Vector vec;
    initVecRandom(&vec, 25);
    Matrix transpose;

    eigenVectorSparse(&mat, &vec, 10);

    printf("Sparse: \n");
    parentArrayPrint(&mat);
    printf("Eigenvector: \n");
    printVec(&vec);

    double eigenValue = rayleighQuotient(&mat, &vec);

    printf("Spectral Radius: %lf\n", eigenValue);

    printf("Verification:\n");

    Vector check;
    initVec(&check, vec.size);
    sparseVecDot(&mat, &vec, &check);
    scalarVec(&vec, eigenValue);

    printf("Altered Eigenvector:\n");
    printVec(&vec);
    printf("Sparse by Eigenvector:\n");
    printVec(&check);*/

    /*eigenVector(&mat, &vec, 10);

    printf("The eigenvector of \n");
    parentArrayPrint(&mat);
    printf("is \n");
    printMat(&vec);*/
    

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
