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

    /*Sparse sparse;
    initSparse(&sparse, 6, 6, 0.35);
    Matrix arnSQ, arnSQT, arnSH;
    initMat(&arnSQ, sparse.rows, sparse.cols);
    initMat(&arnSH, sparse.rows+1, sparse.cols);

    Matrix mat;
    initMat(&mat, sparse.rows, sparse.cols);*/

    /*set(&mat, 2.11, 1, 1);
    set(&mat, -0.27, 3, 1);
    set(&mat, -0.19, 4, 1);

    set(&mat, -0.89, 0, 2);

    set(&mat, -0.52, 1, 3);

    set(&mat, 0.12, 0, 4);
    set(&mat, 0.64, 3, 4);
    set(&mat, -0.71, 4, 4);*/
    
    ESN net;
    initNet(&net, 0, 20, 1, 300, 0.6, 50);
    printWeights(net.weights);
    train(&net);
    //dampening(&net);
    test(&net);

    /*Matrix a, b;
    initMat(&a, sparse.rows, sparse.cols);
    initMat(&b, sparse.rows, sparse.cols);
    
    int o = 1;
    while (arnoldiSparse(&sparse, &arnSQ, &arnSH) != 0) {
        printf("Iteration %d\n", o);
        cleanSparse(&sparse);
        initSparse(&sparse, sparse.rows, sparse.cols, sparse.density);
        cleanMat(&arnSQ);
        initMat(&arnSQ, sparse.rows, sparse.cols);
        cleanMat(&arnSH);
        initMat(&arnSH, sparse.rows+1, sparse.cols);
        o++;
    }

    sparseToMat(&sparse, &mat);
    shrinkMat(&arnSH, arnSQ.rows, arnSQ.cols);
    initTranspose(&arnSQ, &arnSQT);
    
    printf("Sparse:\n");
    printMat(&mat);
    printf("Hessenberg:\n");
    printMat(&arnSH);
    sparseToMat(&sparse, &mat);
    
    matDot(&arnSQT, &mat, &a);
    printf("Check:\n");
    matDot(&a, &arnSQ, &b);
    printMat(&b);
    
    *//*Matrix qr;
    initMat(&qr, 6, 6);
    
    qr.array[0][0] = 7;
    qr.array[0][1] = 7.2761;
    qr.array[0][2] = 5.8120;
    qr.array[0][3] = -0.1397;
    qr.array[0][4] = 9.0152;
    qr.array[0][5] = 7.9363;

    qr.array[1][0] = 12.3693;
    qr.array[1][1] = 4.1307;
    qr.array[1][2] = 18.9685;
    qr.array[1][3] = -1.2071;
    qr.array[1][4] = 10.6833;
    qr.array[1][5] = 2.4160;

    qr.array[2][1] = -7.1603;
    qr.array[2][2] = 2.4478;
    qr.array[2][3] = -0.5656;
    qr.array[2][4] = -4.1814;
    qr.array[2][5] = -3.2510;

    qr.array[3][2] = -8.5988;
    qr.array[3][3] = 2.9151;
    qr.array[3][4] = -3.4169;
    qr.array[3][5] = 5.7230;

    qr.array[4][3] = 1.0464;
    qr.array[4][4] = -2.8351;
    qr.array[4][5] = -10.9792;

    qr.array[5][4] = 1.4143;
    qr.array[5][5] = 5.3415;

    qrHess(&qr);
    printf("QR:\n");
    printMat(&qr);*/

    //printf("QR Convergence:\n");
    //printMat(&qr);

    /*sparseToMat(&sparse, &mat);
    
    
    matDot(&mat, &arnSQ, &a);
    printf("Check:\n");
    matDot(&arnSQ, &arnSH, &b);
    printMat(&b);

    Matrix temp, out, arnD, arnDH, arnDT;
    initMat(&arnD, mat.rows, mat.cols);
    initMat(&arnDH, mat.rows+1, mat.cols);
    initMat(&out, mat.rows, mat.cols);

    arnoldiDense(&mat, &arnD, &arnDH);
    shrinkMat(&arnDH, mat.rows, mat.cols);

    printf("Dense:\n");
    printMat(&mat);
    printf("Hessenberg:\n");
    printMat(&arnDH);

    initTranspose(&arnD, &arnDT);
    initMat(&temp, mat.rows, mat.cols);
    matDot(&arnDT, &mat, &temp);
    matDot(&temp, &arnD, &out);

    printf("Check:\n");
    printMat(&out);*/

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
