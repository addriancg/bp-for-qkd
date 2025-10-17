#include "getchecknodetable.h"

double getRealRate (double rate);
int* matrixToVector(int **matrix, int rows, int cols);

SparseMatrix* generate_P_transpose(double rate);
SparseMatrix* generate_bidiagonal(int m);
SparseMatrix* concatenate_matrices(SparseMatrix *H1, SparseMatrix *H2);

LDPCCode dvbs2ldpcShort(double rate);
SparseMatrix makeSmallTestH();
void freeSparseMatrix(SparseMatrix *H);