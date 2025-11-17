#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "structs.h"
#include "bp_decoder.h"

// Genera n_flips índices distintos en [0, n)
static void generate_random_flips(int n, int n_flips, int *flip_idx, unsigned seed)
{
    int *used = (int*)calloc(n, sizeof(int));
    if (!used) { fprintf(stderr, "OOM in generate_random_flips\n"); exit(1); }
    srand(seed);
    int count = 0;
    while (count < n_flips) {
        int r = rand() % n;
        if (!used[r]) { used[r] = 1; flip_idx[count++] = r; }
    }
    free(used);
}

// ordenar para imprimir bonito
static int cmp_int(const void *a, const void *b){
    int ia = *(const int*)a, ib = *(const int*)b;
    return (ia > ib) - (ia < ib);
}

int main(int argc, char *argv[]){

    unsigned seed = 0; // valor por defecto
    if (argc > 1) {
        seed = (unsigned)atoi(argv[1]);
    }

    // 1) Código DVB-S2 short 1/2
    LDPCCode code = dvbs2ldpcShort(1.0/2.0);

    // 2) Construir decoder (topología y buffers)
    MinSumDecoder *d = create_decoder(&code);
    if (!d || !sanity_check(d)) { fprintf(stderr,"decoder build failed\n"); return 1; }

    // 3) Preparar LLRs con flips aleatorios reproducibles
    double *Lch = (double*)malloc(code.n * sizeof(double));
    if (!Lch) { fprintf(stderr,"OOM Lch\n"); return 1; }

    const double A = 0.8;       // fiabilidad (bajar para más difícil)
    const int n_flips = 30;     // nº de bitflips (subir para más difícil)
    int *flips = (int*)malloc(n_flips * sizeof(int));
    if (!flips) { fprintf(stderr,"OOM flips\n"); return 1; }

    generate_random_flips(code.n, n_flips, flips, seed);
    qsort(flips, n_flips, sizeof(int), cmp_int);

    make_LLRs_static(code.n, A, flips, n_flips, Lch);

    // Info de entrada
    printf(">> Iniciando decoder LDPC (DVB-S2 short)\n");
    printf("   Topología verificada: %s\n", sanity_check(d) ? "OK" : "FAIL");

    printf("\n>> Configuración del canal\n");
    printf("   Fiabilidad (A): %.3f\n", A);
    printf("   Bit flips introducidos: %d\n", n_flips);
    printf("   Semilla: %u\n", seed);
    printf("   Posiciones de flips: ");
    for (int i = 0; i < n_flips; i++) printf("%d%s", flips[i], (i+1<n_flips)?", ":"\n");

    // 4) Métrica del canal (hard decision antes de decodificar)
    int ch_err = 0;
    for (int i = 0; i < code.n; i++) {
        int ch_bit = (Lch[i] < 0.0) ? 1 : 0; // comparar contra todo-ceros
        if (ch_bit) ch_err++;
    }
    printf("\n>> Métricas previas a la decodificación\n");
    printf("   Errores canal (hard decision): %d\n", ch_err);

    // 5) Decodificar: layered normalized min-sum (parámetros “medios”)
    const int max_iter = 30;
    const double alpha = 1.0;   // 1.0 => min-sum normalizado
    const double beta  = 1.0;   // 1.0 => sin damping adicional
    int it = decode(d, &code, Lch, max_iter, alpha, beta);

    // 6) Validación y salida de demo (sin imprimir bits)
    harden(d);
    int ok = check_syndrome(&code, d->hard_decisions);

    int out_err = 0;
    for (int i = 0; i < code.n; i++)
        if (d->hard_decisions[i]) out_err++;

    printf("\n>> Decodificación\n");
    printf("   Parámetros: max_iter=%d, alpha=%.2f, beta=%.2f\n", max_iter, alpha, beta);
    printf("   Iteraciones realizadas: %d\n", (it > 0) ? it : max_iter);

    printf("\n>> Resultados\n");
    printf("   Errores tras decodificación (hard decision): %d\n", out_err);
    printf("   Síndrome final: %s\n", ok ? "OK" : "FAIL");

    int corrected = (ch_err >= out_err) ? (ch_err - out_err) : 0;
    double corr_rate = (ch_err > 0) ? (100.0 * corrected / ch_err) : 100.0;
    printf("   Corrección lograda: %d de %d (%.2f%%)\n", corrected, ch_err, corr_rate);

    if (ok && out_err == 0) {
        printf("\n>> Conclusión: Decodificación completada con éxito. Mensaje recuperado sin errores residuales.\n");
    } else if (ok) {
        printf("\n>> Conclusión: Síndrome OK, quedan %d errores duros. Ajustar fiabilidad o presupuesto de iteraciones.\n", out_err);
    } else {
        printf("\n>> Conclusión: No se alcanzó convergencia en el presupuesto de iteraciones.\n");
    }

    free(flips);
    free(Lch);
    free_decoder(d);
    return ok ? 0 : 2;
}
