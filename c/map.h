/*
 * Author: Joey Tan
 * Date Created: 2-22-20
 * Last Edit: 2-23-20, Joey Tan
 */ 

#ifndef MAP
#define MAP
#include <stdio.h>

// TODO: Renaming sparse matrix representation
// These three structs compose the sparse matrix representation

// The given row and associated value
typedef struct {
    int first;
    double second;
} Pair;

// A column with its filled-in rows
typedef struct {
    Pair* array;
    int value;
    int arraySize;
} Parent;

// High level view of the sparse matrix
// Contains info on the matrix
// as well as the data itself
typedef struct {
    Parent* array;
    int arraySize;
    int cols;
    int rows;
    double density;
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
    pa->arraySize = 0;
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

// These two ensure that added values to the arrays
// are in an ascending sort order
// Uniqueness is handled in initSparse(...)

// add a value to the parent's child array
void parentAdd(Parent* p, int first, double second) {
    if (p->arraySize == 0) {
        p->array = calloc(2, sizeof(Pair));
        p->array[0].first = first;
        p->array[0].second = second;
        p->arraySize++;
    }
    else {
        int index = p->arraySize;
        for (int i = 0; i < p->arraySize; i++) {
            if (first < p->array[i].first) {
                index = i;
                break;
            }
        }

        p->arraySize++;
        Pair* newArray = calloc(p->arraySize, sizeof(Pair));

        for (int i = 0; i < index; i++)
            newArray[i] = p->array[i];
        
        newArray[index].first = first;
        newArray[index].second = second;
        
        for(int i = index+1; i < p->arraySize; i++) {
            newArray[i] = p->array[i-1];
        }

        free(p->array);
        p->array = newArray;
    }
}

void parentArrayAdd(ParentArray* pa, Parent* p) {
    if (pa->arraySize == 1) {
        pa->array = calloc(2, sizeof(Parent));
        pa->array[0] = *p;
        pa->arraySize++;
    }
    else {
        int index = pa->arraySize-1;
        for (int i = 0; i < pa->arraySize; i++) {
            if (p->value < pa->array[i].value) {
                index = i;
                break;
            }
        }

        pa->arraySize++;
        Parent* newArray = calloc(pa->arraySize, sizeof(Parent));

        for (int i = 0; i < index; i++)
            newArray[i] = pa->array[i];
        
        newArray[index] = *p;
        
        for(int i = index+1; i < pa->arraySize; i++) {
            newArray[i] = pa->array[i-1];
        }

        free(pa->array);
        pa->array = newArray;
    }
}

//------------------------------------\\
// I/O                                \\
//------------------------------------\\

void parentPrint(Parent* p) {
    printf("{ Column: %d }\n", p->value);
    for (int i = 0; i < p->arraySize; i++)
        printf("\tRow: %d\n\tValue: %.18lf\n", p->array[i].first, p->array[i].second);
}

void parentArrayPrint(ParentArray* pa) {
    printf("Rows: %d\n", pa->rows);
    printf("Columns: %d\n", pa->cols);
    printf("Density: %lf\n", pa->density);
    printf("[ \n");
    for (int i = 0; i < pa->arraySize-1; i++) {
        printf("\t");
        parentPrint(&(pa->array[i]));
        printf("\n");
    }
    printf("]\n");
}

// TODO: comment on what these actually do

int parentArrayIn(ParentArray* pa, int value) {
    int output = 0;
    for (int i = 0; i < pa->arraySize && pa->array[i].value != value; i++)
        output = i;
    return output;
}

int parentInFirst(Parent* p, int value) {
    int output = 0;
    for (int i = 0; i < p->arraySize; i++)
        if (p->array[i].first == value) return ++output;
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
