#include "structs.h"
#include <stdlib.h>
#include <stdio.h>

int* generate_random_info_bits(int length);
int validate_codeword(const SparseMatrix *H, const int *codeword);
int* encode_ldpc_dvbs2_short(LDPCCode *code, int* info_bits);