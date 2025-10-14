#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structs.h"

typedef struct {
    // Conectividad del grafo (memoria exacta)
    int **var_to_checks;     // var_to_checks[var] = array de checks conectados
    int **check_to_vars;     // check_to_vars[check] = array de vars conectadas
    int *var_degrees;        // grado de cada variable
    int *check_degrees;      // grado de cada check

    int **var_to_checks_pos;
    int **check_to_vars_pos;
    
    // Mensajes (memoria exacta por grado)
    double **v2c_messages;   // v2c_messages[var] = array de mensajes a checks   ###### OPCIONAL ¿?¿?¿ se pueden calcular al vuelo (revisar)
    double **c2v_messages;   // c2v_messages[check] = array de mensajes a vars
    
    // Arrays auxiliares
    double *channel_llrs;    // LLRs del canal [n]
    double *total_llrs;      // LLRs totales [n] 
    int *hard_decisions;     // decisiones hard [n]
    
    // Dimensiones
    int n, m;                // variables y checks
} MinSumDecoder;


void init_frame(MinSumDecoder *d, const double *Lch_in)
{
    // 1) Copiar LLRs del canal
    memcpy(d->channel_llrs, Lch_in, d->n * sizeof(double));

    // 2) Inicializar LLRs totales = canal
    memcpy(d->total_llrs, d->channel_llrs, d->n * sizeof(double));

    // 3) Resetear mensajes check→var a 0 (condición inicial)
    for (int c = 0; c < d->m; c++) {
        int deg = d->check_degrees[c];
        memset(d->c2v_messages[c], 0, deg * sizeof(double));
    }

    // 4) (Opcional) Resetear mensajes var→check si los usas/depuras
    if (d->v2c_messages) {
        for (int v = 0; v < d->n; v++) {
            int deg = d->var_degrees[v];
            memset(d->v2c_messages[v], 0, deg * sizeof(double));
        }
    }

    // 5) (Opcional) Limpiar decisiones duras
    if (d->hard_decisions) {
        memset(d->hard_decisions, 0, d->n * sizeof(int));
    }

    // 6) (Opcional) Rellenar con NaN para detectar usos indebidos durante debug
    // for (int c = 0; c < d->m; c++) {
    //     int deg = d->check_degrees[c];
    //     for (int k = 0; k < deg; k++) d->c2v_messages[c][k] = 0.0; // o NAN
    // }
}

// Construye LLRs estáticos a partir de un patrón de fiabilidad A y un set de posiciones “flip”.
void make_LLRs_static(int n, double A, const int *flip_idx, int n_flips, double *Lch){
    // Todo-ceros: LLR positivo favorece 0
    for (int i = 0; i < n; i++) Lch[i] = +A;
    for (int j = 0; j < n_flips; j++){
        int i = flip_idx[j];
        if (0 <= i && i < n) Lch[i] = -A; // simula bit recibido como 1 con misma fiabilidad
    }
}


MinSumDecoder* create_decoder(LDPCCode *code) {
    MinSumDecoder *decoder = malloc(sizeof(MinSumDecoder));
    if (!decoder) return NULL;
    
    decoder->n = code->n;
    decoder->m = code->m;
    
    // ──── FASE 1: Contar grados ────
    decoder->var_degrees = calloc(decoder->n, sizeof(int));
    decoder->check_degrees = calloc(decoder->m, sizeof(int));
    
    for (int i = 0; i < code->H.num_elements; i++) {
        int check = code->H.row_indices[i];
        int var = code->H.col_indices[i];
        decoder->var_degrees[var]++;
        decoder->check_degrees[check]++;
    }
    
    // ──── FASE 2: Alocar arrays de conectividad (tamaño exacto) ────
    decoder->var_to_checks = malloc(decoder->n * sizeof(int*));
    decoder->check_to_vars = malloc(decoder->m * sizeof(int*));
    decoder->var_to_checks_pos = malloc(decoder->n * sizeof(int*));
    decoder->check_to_vars_pos = malloc(decoder->m * sizeof(int*));
    decoder->v2c_messages = malloc(decoder->n * sizeof(double*));
    decoder->c2v_messages = malloc(decoder->m * sizeof(double*));
    
    // Para cada variable: malloc exacto según su grado
    for (int var = 0; var < decoder->n; var++) {
        int degree = decoder->var_degrees[var];
        decoder->var_to_checks[var] = malloc(degree * sizeof(int));
        decoder->v2c_messages[var] = calloc(degree, sizeof(double));
        decoder->var_to_checks_pos[var] = malloc(degree * sizeof(int));
    }
    
    // Para cada check: malloc exacto según su grado  
    for (int check = 0; check < decoder->m; check++) {
        int degree = decoder->check_degrees[check];
        decoder->check_to_vars[check] = malloc(degree * sizeof(int));
        decoder->c2v_messages[check] = calloc(degree, sizeof(double));
        decoder->check_to_vars_pos[check] = malloc(degree * sizeof(int));
    }
    
    // ──── FASE 3: Construir índices de conectividad ────
    int *var_counters = calloc(decoder->n, sizeof(int));
    int *check_counters = calloc(decoder->m, sizeof(int));
    
    for (int i = 0; i < code->H.num_elements; i++) {
        int check = code->H.row_indices[i];
        int var = code->H.col_indices[i];
        
        int t = var_counters[var]++;
        int k = check_counters[check]++;

        decoder->var_to_checks[var][t] = check;
        decoder->check_to_vars[check][k] = var;
        
        decoder->var_to_checks_pos[var][t] = k;     // en 'check', esta var está en posición k
        decoder->check_to_vars_pos[check][k] = t;   // en 'var', este check está en posición t
        
    }
    
    free(var_counters);
    free(check_counters);
    
    // ──── FASE 4: Arrays auxiliares ────
    decoder->channel_llrs = malloc(decoder->n * sizeof(double));
    decoder->total_llrs = malloc(decoder->n * sizeof(double));
    decoder->hard_decisions = malloc(decoder->n * sizeof(int));
    
    return decoder;
}

void decode_one_iteration_layered(MinSumDecoder *d, double alpha, double beta)
{
    // alpha: factor de normalized min-sum (1.0 => min-sum clásico)
    // beta : damping para C2V: new = (1-beta)*old + beta*raw  (1.0 => sin damping)

    for (int c = 0; c < d->m; c++) {
        const int deg = d->check_degrees[c];
        int *nbrs = d->check_to_vars[c];

        // Buffers temporales por check (evita malloc para grados pequeños)
        double v2c_vals_stack[64];
        double abs_stack[64];
        int    sgn_stack[64];

        double *v2c_vals = (deg <= 64) ? v2c_vals_stack : (double*)malloc(deg*sizeof(double));
        double *abs_vals = (deg <= 64) ? abs_stack     : (double*)malloc(deg*sizeof(double));
        int    *sgn_vals = (deg <= 64) ? sgn_stack     : (int*)malloc(deg*sizeof(int));

        int sign_global = +1;
        double min1 = INFINITY, min2 = INFINITY;
        int idx_min1 = -1;

        // Primer barrido: construir v2c, signo global y min1/min2
        for (int k = 0; k < deg; k++) {
            int v = nbrs[k];
            double c2v_old = d->c2v_messages[c][k];
            double v2c = d->total_llrs[v] - c2v_old;    // extrínseco al vuelo

            v2c_vals[k] = v2c;
            double a = fabs(v2c);
            abs_vals[k] = a;

            int s = (v2c >= 0.0) ? +1 : -1;
            sgn_vals[k] = s;
            sign_global *= s;

            if (a < min1) {
                min2 = min1;
                min1 = a;
                idx_min1 = k;
            } else if (a < min2) {
                min2 = a;
            }
        }

        // Segundo barrido: formar c2v nuevos y actualizar total_llrs in-place
        for (int k = 0; k < deg; k++) {
            int v = nbrs[k];

            double mag = (k == idx_min1) ? min2 : min1;
            int sign_out = sign_global * sgn_vals[k];    // excluye al propio k

            double raw = alpha * (double)sign_out * mag; // alpha=1.0 => min-sum
            double old = d->c2v_messages[c][k];
            double newv = (1.0 - beta) * old + beta * raw; // beta=1.0 => sin damping

            d->total_llrs[v] += (newv - old);
            d->c2v_messages[c][k] = newv;

            if (d->v2c_messages) {
                int t = d->check_to_vars_pos[c][k];
                d->v2c_messages[v][t] = v2c_vals[k];
            }
        }

        if (v2c_vals != v2c_vals_stack) free(v2c_vals);
        if (abs_vals != abs_stack) free(abs_vals);
        if (sgn_vals != sgn_stack) free(sgn_vals);
    }
}

void harden(const MinSumDecoder *d)
{
    for (int v = 0; v < d->n; v++) {
        d->hard_decisions[v] = (d->total_llrs[v] < 0.0) ? 1 : 0;
    }
}

// Asumimos H en formato lista de 1s: code->H.num_elements pares (row,col)
int check_syndrome(const LDPCCode *code, const int *hard_bits)
{
    int m = code->m;
    int ok = 1;
    int *acc = (int*)calloc(m, sizeof(int));
    for (int i = 0; i < code->H.num_elements; i++) {
        int r = code->H.row_indices[i];
        int c = code->H.col_indices[i];
        acc[r] ^= hard_bits[c]; // XOR
    }
    for (int r = 0; r < m; r++) {
        if (acc[r] & 1) { ok = 0; break; }
    }
    free(acc);
    return ok;
}

int decode(MinSumDecoder *d, const LDPCCode *code,
           const double *Lch_in, int max_iter,
           double alpha, double beta)
{
    init_frame(d, Lch_in);

    for (int it = 1; it <= max_iter; ++it) {
        decode_one_iteration_layered(d, alpha, beta);
        harden(d);
        if (check_syndrome(code, d->hard_decisions)) {
            return it; // éxito
        }
    }
    return -1; // no convergió en max_iter
}

int sanity_check(const MinSumDecoder *d){
    for (int v = 0; v < d->n; v++){
        for (int t = 0; t < d->var_degrees[v]; t++){
            int c = d->var_to_checks[v][t];
            int k = d->var_to_checks_pos[v][t];
            if (d->check_to_vars[c][k] != v) return 0;
            if (d->check_to_vars_pos[c][k] != t) return 0;
        }
    }
    return 1;
}

void free_decoder(MinSumDecoder *d){
    if (!d) return;
    for (int v = 0; v < d->n; v++){
        free(d->var_to_checks[v]);
        free(d->var_to_checks_pos[v]);
        free(d->v2c_messages[v]);
    }
    for (int c = 0; c < d->m; c++){
        free(d->check_to_vars[c]);
        free(d->check_to_vars_pos[c]);
        free(d->c2v_messages[c]);
    }
    free(d->var_to_checks); free(d->var_to_checks_pos);
    free(d->check_to_vars); free(d->check_to_vars_pos);
    free(d->v2c_messages);  free(d->c2v_messages);
    free(d->channel_llrs);  free(d->total_llrs);
    free(d->hard_decisions);
    free(d);
}
