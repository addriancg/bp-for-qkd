#pragma once

typedef struct MatrixInt MatrixInt;
typedef struct Matrices Matrices;
typedef struct SparseMatrix SparseMatrix;
typedef struct SparseMatrixCSR SparseMatrixCSR;
typedef struct LDPCCode LDPCCode;

struct MatrixInt
{
	int **data;
    int rows;
    int cols;
};

struct Matrices
{
	MatrixInt ct1;
	MatrixInt ct2;
};

struct SparseMatrix
{
	int *row_indices;
	int *col_indices;
	int *values;
	int num_elements;
	int num_rows;
	int num_cols;
};

struct SparseMatrixCSR
{
    int *indptr;
    int *indices;
    int *values;
    int num_elements;
	int num_rows;
	int num_cols;
};

struct LDPCCode
{
    SparseMatrix H;
    int k;
    int n;
    int m;
};
