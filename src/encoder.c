#include "encoder.h"


/* 	Codifica info_bits[] de longitud k, devuelve codeword[] de longitud n (n = k + m)
	Teniendo en cuenta que se está siguiendo el estandar  DVB-S2 short la estructura es
	sistematica: codeword = [ info  ||  paridad ]
*/

int* generate_random_info_bits(int length){
	int *bits = malloc(length * sizeof(int));
    if (!bits) {
        fprintf(stderr, "Error: malloc failed in generate_random_bits\n");
        return NULL;
    }
    
    for (int i = 0; i < length; i++) {
        bits[i] = rand() & 1;
    }

	return bits;
	
}

int validate_codeword(const SparseMatrix *H, const int *codeword) {
    for (int row = 0; row < H->num_rows; row++) {
        int syndrome = 0;
        for (int idx = H->row_indices[row]; idx < H->row_indices[row + 1]; idx++) {
            syndrome ^= codeword[H->col_indices[idx]];
        }
        if (syndrome != 0) return 0;
    }
    return 1;
}

int* encode_ldpc_dvbs2_short(SparseMatrix *H, int* info_bits, int k, int m){
	int n = k + m;
    int *codeword = malloc(n * sizeof(int));
	if (!codeword) {
        fprintf(stderr, "Error: malloc failed in encode_ldpc_dvbs2_short\n");
        return NULL;
    }

	for (int i = 0; i < k; i++) {
        codeword[i] = info_bits[i];
    }

	for (int row = 0; row < m; row++) {
        int parity = 0;
		
		for (int idx = H->row_indices[row]; idx < H->row_indices[row + 1]; idx++) {
            int col = H->col_indices[idx];
            if (col < k) {
                parity ^= codeword[col]; // XOR de los bits de información
            }
        }

		if (row > 0) {
            parity ^= codeword[k + row - 1];
        }

        codeword[k + row] = parity;
    }

	return codeword;
}