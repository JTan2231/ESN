#include <stdio.h>
#include <stdlib.h>

typedef struct PairTag {
    int first;
    double second;
} Pair;

typedef struct ColumnTag {
    Pair* array;
    int value;
    int arraySize;
} Column;

typedef struct SparseTag {
    Column* array;
    int arraySize;
    int cols;
    int rows;
    double density;
} Sparse;

void initColumn(Column* p, int value) {
    p->array = malloc(sizeof p->array);
    p->array[0].first = -1;
    p->value = value;
    p->arraySize = 0;
}

void cleanColumn(Column* p) {
    free(p->array);
    p->arraySize = 0;
}

void cleanSparse(Sparse* pa) {
    for (int i = 0; i < pa->arraySize; i++)
        cleanColumn(&(pa->array[i]));
    free(pa->array);
    pa->arraySize = 0;
}

// These two ensure that added values to the arrays
// are in an ascending sort order
// Uniqueness is handled in initSparse(...)

// add a value to the Column's child array
void columnAdd(Column* p, int first, double second) {
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
        Pair* newArray = calloc(p->arraySize+1, sizeof(Pair));

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

// add new column to a sparse matrix
void sparseAdd(Sparse* pa, Column* p) {
    if (pa->arraySize == 1) {
        pa->array = calloc(2, sizeof(Column));
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
        Column* newArray = calloc(pa->arraySize+1, sizeof(Column));

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

void columnPrint(Column* p) {
    printf("{ Column: %d }\n", p->value);
    for (int i = 0; i < p->arraySize; i++)
        printf("\tRow: %d\n\tValue: %.18lf\n", p->array[i].first, p->array[i].second);
}

void sparsePrint(Sparse* pa) {
    printf("Rows: %d\n", pa->rows);
    printf("Columns: %d\n", pa->cols);
    printf("Density: %lf\n", pa->density);
    printf("Units: %d\n", (int)(pa->rows*pa->cols*pa->density));
    printf("[ \n");
    for (int i = 0; i < pa->arraySize-1; i++) {
        printf("\t");
        columnPrint(&(pa->array[i]));
        printf("\n");
    }
    printf("]\n");
}

// TODO: comment on what these actually do

int sparseIn(Sparse* pa, int value) {
    int output = 0;
    for (int i = 0; i < pa->arraySize && pa->array[i].value != value; i++)
        output = i;
    return output;
}

int columnInFirst(Column* p, int value) {
    int output = 0;
    for (int i = 0; i < p->arraySize; i++)
        if (p->array[i].first == value) return ++output;
    return output;
}

int columnInSecond(Column* p, int value) {
    int output = 0;
    for (int i = 0; i < p->arraySize && p->array[i].second != value; i++)
        output = i;
    return output;
}

int columnIn(Column* p, int first, int second) {
    int output = 0;
    for (int i = 0;
         i < p->arraySize && p->array[i].first != first && p->array[i].second != second;
         i++)
        output = i;
    return output;
}
