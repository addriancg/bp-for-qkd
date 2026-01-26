#ifndef BP_DECODER_H
#define BP_DECODER_H


#include "matrixUtils.h"
#include "dvbs2ldpcShort.h"
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

    //Sindrome
    int *syndrome_received;
} MinSumDecoder;




// Funciones principales
MinSumDecoder* create_decoder(LDPCCode *code);
void free_decoder(MinSumDecoder *d);
void make_LLRs_static(int n, double A, const int *flip_idx, int n_flips, double *Lch);
int decode(MinSumDecoder *d, const LDPCCode *code, const double *Lch_in, int max_iter, double alpha, double beta);
void harden(const MinSumDecoder *d);
int check_syndrome(const LDPCCode *code, const int *hard_bits);
int check_syndrome_match(const LDPCCode *code, const int *hard_bits, const int *syndrome_received);
int sanity_check(const MinSumDecoder *d);

#endif /* BP_DECODER_H */