/*
 * Author: Joey Tan
 * Date Created: 2-22-20
 * Last Edit: 2-22-20, Joey Tan
 */ 

#include <stdio.h>

typedef struct {
    int* array;
    int value;
    int arraySize;
} Parent;

typedef struct {
    Parent* array;
    int arraySize;
} ParentArray;

//------------------------------------\\
// Initialization                     \\
//------------------------------------\\

void initParent(Parent* p, int value) {
    p->array = malloc(sizeof p->array);
    p->array[0] = -1;
    p->value = value;
    p->arraySize = 1;
}

void initParentArray(ParentArray* pa) {
    pa->array = malloc(sizeof pa->array);
    pa->arraySize = 1;
}

//------------------------------------\\
// Cleanup                            \\
//------------------------------------\\

void cleanParent(Parent* p) {
    free(p->array);
    p->arraySize = 0;
}

void cleanParentArray(ParentArray* pa) {
    for (int i = 0; i < pa->arraySize; i++)
        cleanParent(&(pa->array[i]));
    free(pa->array);
    pa->arraySize = 0;
}

//------------------------------------\\
// Mutation                           \\
//------------------------------------\\

// add a value to the parent's child array
void parentAdd(Parent* p, int value) {
    p->array[p->arraySize-1] = value;
    p->arraySize++;

    // expand array
    int* newArray = malloc(sizeof p->array * p->arraySize);
    for (int i = 0; i < p->arraySize; i++)
        newArray[i] = p->array[i];
    free(p->array);
    p->array = newArray;
    p->array[p->arraySize] = -1;
}

void parentArrayAdd(ParentArray* pa, Parent* p) {
    pa->array[pa->arraySize-1] = *p;
    pa->arraySize++;

    // expand array
    Parent* newArray = malloc(sizeof pa->array * pa->arraySize);
    for (int i = 0; i < pa->arraySize; i++)
        newArray[i] = pa->array[i];
    free(pa->array);
    pa->array = newArray;
}

//------------------------------------\\
// I/O                                \\
//------------------------------------\\

void parentPrint(Parent* p) {
    printf("[ ");
    for (int i = 0; i < p->arraySize; i++)
        printf("%d ", p->array[i]);
    printf("]\n");
}

void parentArrayPrint(ParentArray* pa) {
    printf("[ \n");
    for (int i = 0; i < pa->arraySize; i++) {
        printf("    ");
        parentPrint(&(pa->array[i]));
        printf("\n");
    }
    printf("]\n");
}
