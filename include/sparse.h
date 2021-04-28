#ifndef MAP
#define MAP
#include <stdio.h>

// TODO: Renaming sparse matrix representation
// These three structs compose the sparse matrix representation

// The given row and associated value
typedef struct PairTag {
    int first;
    double second;
} Pair;

// A column with its filled-in rows
// Column.array is a dynamic array
// with key-value entries (row-value)
// added as the matrix is constructed
//
// rows for which there is no key (e.g. row[i])
// are interpreted as being zero
// at point (row[i], Column.value)
typedef struct ColumnTag {
    Pair* array;
    int value;
    int arraySize;
} Column;

// Sparse matrix representation
//
// Memory is only used on matrix info
// and points at which there is a
// non-zero number
//
// Contains an array of columns
// in which are key-value pairs
// indicating where the non-zero
// values are (i.e. which row
//                  in that column)
//
// density == (# of non-zero values) / (rows * cols)
typedef struct SparseTag {
    Column* array;
    int arraySize;
    int cols;
    int rows;
    double density;
} Sparse;

//------------------------------------\\
// Initialization                     \\
//------------------------------------\\

void initColumn(Column* p, int value);

//------------------------------------\\
// Cleanup                            \\
//------------------------------------\\

void cleanColumn(Column* p);

void cleanSparse(Sparse* pa);

//------------------------------------\\
// Mutation                           \\
//------------------------------------\\

// These two ensure that added values to the arrays
// are in an ascending sort order
// Uniqueness is handled in initSparse(...)

// add a value to the Column's child array
void columnAdd(Column* p, int first, double second);

// add new column to a sparse matrix
void sparseAdd(Sparse* pa, Column* p);

//------------------------------------\\
// I/O                                \\
//------------------------------------\\

void columnPrint(Column* p);

void sparsePrint(Sparse* pa);

// TODO: comment on what these actually do

int sparseIn(Sparse* pa, int value);

int columnInFirst(Column* p, int value);

int columnInSecond(Column* p, int value);

int columnIn(Column* p, int first, int second);

#endif
