#pragma once
#include <stdint.h>
#include "structs.h"

int* addcr(const int *c, const int *r, int M, int N);
int* matrixToVector(int **matrix, int rows, int cols);
void freeMatrix(int **matrix, int rows);
void printMatrix(int **matriz, int rows, int cols);
SparseMatrix coo_to_csr(int num_rows, int num_cols, int nnz, const int *coo_row, const int *coo_col);
int isSystematicForm(SparseMatrix H, int k);
//SparseLIL csr_to_lil(const SparseMatrix *csr);
//void free_lil(SparseLIL *lil);
//void xor_rows(SparseLIL *lil, int dest_row, int src_row);
//int gaussian_elimination_systematic(SparseLIL *lil, int k, int *col_perm);
//SparseMatrix lil_to_csr(const SparseLIL *lil);
uint8_t **csr_to_dense(const SparseMatrix *csr);
void free_dense(uint8_t **matrix, int rows);
uint8_t **compute_nullspace_gf2(uint8_t **H, int m, int n, int *k_out);
int check_HGt(uint8_t **H, int m, int n, uint8_t **G, int k);
int gcd(int a, int b);
uint8_t **loadG(const char *filename, int *k_out, int *n_out);
int saveG(const char *filename, uint8_t **G, int k, int n);