/*
 * Author: Joey Tan
 * Date Created: 2-22-20
 * Last Edit: 2-23-20, Joey Tan
 */ 

#ifndef MAP
#define MAP
#include <stdio.h>

typedef struct {
    int first;
    float second;
} Pair;

typedef struct {
    Pair* array;
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
    p->array[0].first = -1;
    p->value = value;
    p->arraySize = 0;
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
void parentAdd(Parent* p, int first, float second) {
    p->array[p->arraySize].first = first;
    p->array[p->arraySize].second = second;
    p->arraySize++;

    // expand array
    Pair* newArray = malloc(sizeof(Pair) * p->arraySize);
    for (int i = 0; i < p->arraySize; i++)
        newArray[i] = p->array[i];
    free(p->array);
    p->array = newArray;
    //p->array[p->arraySize].first = 0;
}

void parentArrayAdd(ParentArray* pa, Parent* p) {
    pa->array[pa->arraySize-1] = *p;
    pa->arraySize++;

    // expand array
    Parent* newArray = malloc(sizeof(Parent) * pa->arraySize);
    for (int i = 0; i < pa->arraySize; i++)
        newArray[i] = pa->array[i];
    free(pa->array);
    pa->array = newArray;
}

//------------------------------------\\
// I/O                                \\
//------------------------------------\\

void parentPrint(Parent* p) {
    printf("[ {value: %d} ", p->value);
    for (int i = 0; i < p->arraySize; i++)
        printf("First, Second: %d, %f ", p->array[i].first, p->array[i].second);
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

int parentArrayIn(ParentArray* pa, int value) {
    int output = 0;
    for (int i = 0; i < pa->arraySize && pa->array[i].value != value; i++)
        output = i;
    return output;
}

int parentInFirst(Parent* p, int value) {
    int output = 0;
    for (int i = 0; i < p->arraySize && p->array[i].first != value; i++)
        output = i;
    return output;
}

int parentInSecond(Parent* p, int value) {
    int output = 0;
    for (int i = 0; i < p->arraySize && p->array[i].second != value; i++)
        output = i;
    return output;
}

int parentIn(Parent* p, int first, int second) {
    int output = 0;
    for (int i = 0;
         i < p->arraySize && p->array[i].first != first && p->array[i].second != second;
         i++)
        output = i;
    return output;
}

#endif
