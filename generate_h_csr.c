#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include "bp_decoder.h"
#include "structs.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        errx(EXIT_FAILURE, "Usage: ./generate_h_ldpc <rate>");
    }

    double rate;
    if (strcmp(argv[1], "1/2") == 0) rate = 1.0/2.0;
    else if (strcmp(argv[1], "2/3") == 0) rate = 2.0/3.0;
    else if (strcmp(argv[1], "3/5") == 0) rate = 3.0/5.0;
    else if (strcmp(argv[1], "5/6") == 0) rate = 5.0/6.0;
    else if (strcmp(argv[1], "8/9") == 0) rate = 8.0/9.0;
    else if (strcmp(argv[1], "1/4") == 0) rate = 1.0/4.0;
    else if (strcmp(argv[1], "1/3") == 0) rate = 1.0/3.0;
    else if (strcmp(argv[1], "2/5") == 0) rate = 2.0/5.0;
    else if (strcmp(argv[1], "3/4") == 0) rate = 3.0/4.0;
    else if (strcmp(argv[1], "4/5") == 0) rate = 4.0/5.0;
    else {
        char *endptr;
        rate = strtod(argv[1], &endptr);
        if (*endptr != '\0') {
            fprintf(stderr, "Error: usa '1/2', '2/3', etc. o número decimal\n");
            return 1;
        }
    }

    printf("Rate: %s → %.6f\n", argv[1], rate);

    LDPCCode code = dvbs2ldpcShort(rate);
    SparseMatrixCSR *H = coo_to_csr(&code.H);
    
    if (!H) {
        fprintf(stderr, "Error: coo_to_csr falló\n");
        return 1;
    }

    FILE *f = fopen("H_matrix_csr.txt", "w");
    if (!f) {
        perror("Error abriendo archivo");
        return 1;
    }

    for(int row = 0; row < H->num_rows; row++) {
        int start = H->indptr[row];
        int end = H->indptr[row + 1];
        for(int i = start; i < end; i++) {
            fprintf(f, "%d %d %d\n", row, H->indices[i], H->values[i]);
        }
    }
    
    printf("-> Generated H_matrix_csr.txt (nnz=%d)\n", H->num_elements);
    fclose(f);
    return 0;
}
