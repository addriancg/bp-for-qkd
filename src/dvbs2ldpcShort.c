// dvbs2ldpcShort.c  –  DVB-S2 Short LDPC (versión estable)
// --------------------------------------------------------
//  • Genera H (9000 × 16200) completamente en dominio QC.
//  • Para cada fila de paridad v crea los 2 ‘1’s de la
//    escalera P^d en la MISMA fila pero con shifts distintos:
//         • col base+⌊(v-1)/Z⌋, shift = (v-1) mod Z
//         • col base+⌊ v   /Z⌋, shift =  v    mod Z
//  • Filtra duplicados pares (GF(2)) antes de pasar a CSR.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dvbs2ldpcShort.h"
#include "matrixUtils.h"

// ───── utilidades ─────────────────────────────────────────────────────────
typedef struct { int r, c; } Pair;
static int cmp_pair(const void *a,const void *b){
    const Pair *pa=a,*pb=b;
    return (pa->r!=pb->r) ? pa->r-pb->r : pa->c-pb->c;
}

// ───── tablas de tasa ------------------------------------------------------
static const double rates[]     ={1.0/4,1.0/3,2.0/5,1.0/2,3.0/5,2.0/3,
                                  3.0/4,4.0/5,5.0/6,8.0/9,9.0/10};
static const double realRates[] ={1.0/5,1.0/3,2.0/5,4.0/9,3.0/5,2.0/3,
                                 11.0/15,7.0/9,37.0/45,8.0/9};
 double getReal(double r){
    for(int i=0;i<11;i++) if(r==rates[i]) return realRates[i];
    return -1.0;
}

// ───── memoria thin-wrapper ----------------------------------------------
void freeSparseMatrix(SparseMatrix *H){
    free(H->row_indices); free(H->col_indices); free(H->values);
    memset(H,0,sizeof *H);
}

// ───── constructor principal ──────────────────────────────────────────────
// ═══════════════════════════════════════════════════════════════
// FUNCIÓN 1: Generar Submatriz P^T (Parte Sistemática)
// ═══════════════════════════════════════════════════════════════
SparseMatrix* generate_P_transpose(double rate) {
    const int n = 16200, Z = 360;
    double Rr = getReal(rate);
    int k = (int)(n * Rr), m = n - k;
    int q = m / Z;
    
    // Obtener tablas del estándar
    Matrices M = getchecknodetable(rate);
    
    // Arrays temporales para coordenadas
    int *temp_rows = malloc(50000 * sizeof(int));  // Solo H₁
    int *temp_cols = malloc(50000 * sizeof(int));
    int count = 0;
    
    // Procesar ct1
    int ct1_offset = 0;
    for(int row = 0; row < M.ct1.rows; row++) {
        for(int col = 0; col < M.ct1.cols; col++) {
            int x = M.ct1.data[row][col];
            for(int m_idx = 0; m_idx < Z; m_idx++) {
                int parity_addr = (x + m_idx * q) % m;
                int info_addr = ct1_offset + row * Z + m_idx;
                if(info_addr < k) {
                    temp_rows[count] = parity_addr;
                    temp_cols[count] = info_addr;
                    count++;
                }
            }
        }
    }
    
    // Procesar ct2
    int ct2_offset = M.ct1.rows * Z;
    for(int row = 0; row < M.ct2.rows; row++) {
        for(int col = 0; col < M.ct2.cols; col++) {
            int x = M.ct2.data[row][col];
            for(int m_idx = 0; m_idx < Z; m_idx++) {
                int parity_addr = (x + m_idx * q) % m;
                int info_addr = ct2_offset + row * Z + m_idx;
                if(info_addr < k) {
                    temp_rows[count] = parity_addr;
                    temp_cols[count] = info_addr;
                    count++;
                }
            }
        }
    }
    
    // Crear matriz sparse P^T
    SparseMatrix *P_t = malloc(sizeof(SparseMatrix));
    P_t->num_rows = m;
    P_t->num_cols = k;
    P_t->num_elements = count;
    P_t->row_indices = malloc(count * sizeof(int));
    P_t->col_indices = malloc(count * sizeof(int));
    P_t->values = malloc(count * sizeof(int));
    
    for(int i = 0; i < count; i++) {
        P_t->row_indices[i] = temp_rows[i];
        P_t->col_indices[i] = temp_cols[i];
        P_t->values[i] = 1;
    }
    
    free_matrices(&M);
    free(temp_rows);
    free(temp_cols);
    return P_t;
}

// ═══════════════════════════════════════════════════════════════
// FUNCIÓN 2: Generar Submatriz Bi-diagonal (Parte de Paridad I)
// ═══════════════════════════════════════════════════════════════
SparseMatrix* generate_bidiagonal(int m) {
    // Calcular número de elementos: diagonal + subdiagonal
    int num_elements = m + (m - 1);
    
    // Crear matriz bi-diagonal
    SparseMatrix *H2 = malloc(sizeof(SparseMatrix));
    H2->num_rows = m;
    H2->num_cols = m;
    H2->num_elements = num_elements;
    H2->row_indices = malloc(num_elements * sizeof(int));
    H2->col_indices = malloc(num_elements * sizeof(int));
    H2->values = malloc(num_elements * sizeof(int));
    
    int count = 0;
    
    // Diagonal principal: I[i,i] = 1
    for(int i = 0; i < m; i++) {
        H2->row_indices[count] = i;
        H2->col_indices[count] = i;
        H2->values[count] = 1;
        count++;
    }
    
    // Diagonal inferior: I[i+1,i] = 1
    for(int i = 0; i < m - 1; i++) {
        H2->row_indices[count] = i + 1;
        H2->col_indices[count] = i;
        H2->values[count] = 1;
        count++;
    }
    
    return H2;
}

// ═══════════════════════════════════════════════════════════════
// FUNCIÓN 3: Concatenar Matrices Horizontalmente H = [P | I]
// ═══════════════════════════════════════════════════════════════
SparseMatrix* concatenate_matrices(SparseMatrix *H1, SparseMatrix *H2) {
    // Verificar compatibilidad
    if(H1->num_rows != H2->num_rows) {
        fprintf(stderr, "Error: matrices deben tener el mismo número de filas\n");
        return NULL;
    }
    
    int total_elements = H1->num_elements + H2->num_elements;
    int total_cols = H1->num_cols + H2->num_cols;
    
    // Crear matriz concatenada
    SparseMatrix *H = malloc(sizeof(SparseMatrix));
    H->num_rows = H1->num_rows;
    H->num_cols = total_cols;
    H->num_elements = total_elements;
    H->row_indices = malloc(total_elements * sizeof(int));
    H->col_indices = malloc(total_elements * sizeof(int));
    H->values = malloc(total_elements * sizeof(int));
    
    int count = 0;
    
    // Copiar H₁ (sin offset de columnas)
    for(int i = 0; i < H1->num_elements; i++) {
        H->row_indices[count] = H1->row_indices[i];
        H->col_indices[count] = H1->col_indices[i];
        H->values[count] = H1->values[i];
        count++;
    }
    
    // Copiar H₂ (con offset de columnas)
    for(int i = 0; i < H2->num_elements; i++) {
        H->row_indices[count] = H2->row_indices[i];
        H->col_indices[count] = H2->col_indices[i] + H1->num_cols;  // Offset!
        H->values[count] = H2->values[i];
        count++;
    }
    
    return H;
}

// ═══════════════════════════════════════════════════════════════
// Creacion matriz H
// ═══════════════════════════════════════════════════════════════
LDPCCode dvbs2ldpcShort(double rate) {
    const int n = 16200;
    double Rr = getReal(rate);
    if(Rr < 0) {
        fprintf(stderr, "Error: rate no soportado\n");
        exit(1);
    }
    
    int k = (int)(n * Rr);
    int m = n - k;
    
    // 1) Generar submatriz P^T (parte sistemática)
    SparseMatrix *P_t = generate_P_transpose(rate);
    if(!P_t) {
        fprintf(stderr, "Error generando P^T\n");
        exit(1);
    }
    
    // 2) Generar submatriz bi-diagonal
    SparseMatrix *H2 = generate_bidiagonal(m);
    if(!H2) {
        fprintf(stderr, "Error generando matriz bi-diagonal\n");
        freeSparseMatrix(P_t); free(P_t);
        exit(1);
    }
    
    // 3) Concatenar matrices H = [P^T | H₂]
    SparseMatrix *H = concatenate_matrices(P_t, H2);
    if(!H) {
        fprintf(stderr, "Error concatenando matrices\n");
        freeSparseMatrix(P_t); free(P_t);
        freeSparseMatrix(H2); free(H2);
        exit(1);
    }
    
    // 4) Crear estructura LDPCCode
    LDPCCode result = {
        .H = *H,
        .k = k,
        .n = n,
        .m = m
    };
    
    // 5) Limpiar matrices temporales
    freeSparseMatrix(P_t); free(P_t);
    freeSparseMatrix(H2); free(H2);
    free(H);  // Solo el struct, no el contenido (se movió a result.H)
    
    return result;
}