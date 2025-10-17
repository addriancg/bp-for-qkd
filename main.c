#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "bp_decoder.h"


int main(void){
    // 1) Código DVB-S2 short 1/2
    LDPCCode code = dvbs2ldpcShort(1.0/2.0);

    // 2) Construir decoder (topología y buffers)
    MinSumDecoder *d = create_decoder(&code);
    if (!d || !sanity_check(d)) { fprintf(stderr,"decoder build failed\n"); return 1; }

    // 3) LLRs estáticos reproducibles
    double *Lch = (double*)malloc(code.n * sizeof(double));
    int flips[] = {5, 37, 102, 511, 2047, 4095};
    make_LLRs_static(code.n, 2.5, flips, (int)(sizeof(flips)/sizeof(flips[0])), Lch);

    // 4) Decodificar: layered normalized min-sum
    int it = decode(d, &code, Lch, /*max_iter*/20, /*alpha*/0.9, /*beta*/1.0);

    // 5) Validación y salida mínima
    harden(d);
    int ok = check_syndrome(&code, d->hard_decisions);
    printf("resultado=%s, iter=%d\n", ok? "OK":"FAIL", (it>0)? it:20);
    for (int i=0; i<64 && i<code.n; i++) putchar('0' + d->hard_decisions[i]);
    putchar('\n');

    free(Lch);
    free_decoder(d);
    return ok? 0 : 2;
}
