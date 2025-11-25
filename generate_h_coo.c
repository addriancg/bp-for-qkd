#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include "bp_decoder.h"
#include "structs.h"



int main (int argc, char *argv[])
{

    if (argc < 2) {
        errx(EXIT_FAILURE,"Usage: ./generate_h_ldpc <rate>");
    }

    char *endptr;
    double rate = strtod(argv[1], &endptr);
    if (*endptr != '\0') {
        fprintf(stderr, "Err: argument not a valid number: %s\n", argv[1]);
        return 1;
    }

    printf("Rate: %f\n", rate);

    LDPCCode code = dvbs2ldpcShort(rate);

    FILE *f = fopen("H_matrix_coo.txt", "w");
    if(!f) {
        perror("Error abriendo archivo");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < code.H.num_elements; i++) {
        fprintf(f, "%d %d %d\n", code.H.row_indices[i], code.H.col_indices[i], code.H.values[i]);
    }
    
    printf("-> Generated .txt file\n");

    fclose(f);

}