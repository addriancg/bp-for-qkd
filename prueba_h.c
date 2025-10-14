#include "dvbs2ldpcShort.h"
#include "matrixUtils.h"
#include "bp_decoder.h"
#include "crc.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751
#endif


int main() {
    printf(">>> INICIO prueba_h\n"); fflush(stdout);

    srand(time(NULL));
    printf("1) srand listo\n"); fflush(stdout);

    // 1) Genera H (ya calcula todo internamente)
    printf("=== GENERACIÓN DE MATRIZ H ===\n");
    LDPCCode code = dvbs2ldpcShort(1.0/2.0, "sparse");
    printf("2) H generado: m=%d, k=%d, n=%d\n",
           code.m, code.k, code.n);
    fflush(stdout);

    SparseMatrix H = code.H;

    // 2) Verificaciones usando los valores YA calculados
    printf("\n=== VERIFICACIONES DE LA MATRIZ H ===\n");
    
    // 2.1) Verificar dimensiones básicas
    printf("2.1) Dimensiones H: %dx%d\n", H.num_rows, H.num_cols);
    printf("     Elementos no-zero: %d\n", H.num_elements);
    
    // 2.2) Verificar coherencia k+m=n
    printf("2.2) Coherencia: k=%d + m=%d = %d %s n=%d\n", 
           code.k, code.m, code.k+code.m,
           (code.k+code.m==code.n)?"==":"!=", code.n);
    
    // 2.3) Calcular densidad
    double density = (double)H.num_elements / ((double)H.num_rows * H.num_cols);
    printf("2.3) Densidad: %.4f%% %s\n", 
           density * 100, (density < 0.01) ? "✓ (buena para LDPC)" : "⚠️ alta");
    
    // 2.4) Verificar que no hay índices fuera de rango
    int invalid_indices = 0;
    for(int i = 0; i < H.num_elements; i++) {
        if(H.row_indices[i] >= H.num_rows || H.row_indices[i] < 0 ||
           H.col_indices[i] >= H.num_cols || H.col_indices[i] < 0) {
            invalid_indices++;
        }
    }
    printf("2.4) Índices inválidos: %d %s\n", 
           invalid_indices, (invalid_indices==0)?"✓":"❌");
    
    // 2.5) Verificar separación H₁ y H₂
    int h1_elements = 0, h2_elements = 0;
    for(int i = 0; i < H.num_elements; i++) {
        if(H.col_indices[i] < code.k) {
            h1_elements++;  // Parte sistemática
        } else {
            h2_elements++;  // Parte de paridad
        }
    }
    printf("2.5) Distribución: H₁=%d elementos, H₂=%d elementos\n", 
           h1_elements, h2_elements);
    
    // 2.6) Verificación específica de H₂ (estructura bi-diagonal)
    int h2_expected = code.m + (code.m - 1);  // diagonal + subdiagonal
    printf("2.6) H₂: encontrados=%d, esperados=%d %s\n", 
           h2_elements, h2_expected, 
           (h2_elements == h2_expected) ? "✓" : "⚠️");
    
    // 2.7) Mostrar algunos elementos para debug
    printf("2.7) DEBUG primeros 10 elementos:\n");
    for(int i = 0; i < 10 && i < H.num_elements; i++) {
        printf("     [%d]: fila=%d, col=%d, val=%d\n", 
               i, H.row_indices[i], H.col_indices[i], H.values[i]);
    }
    
    fflush(stdout);

    // 3) Resumen final
    printf("\n=== RESUMEN FINAL ===\n");
    if(code.k + code.m == code.n && invalid_indices == 0 && density < 0.01) {
        printf("✅ MATRIZ H GENERADA CORRECTAMENTE\n");
        printf("   - Dimensiones correctas: %dx%d\n", H.num_rows, H.num_cols);
        printf("   - Densidad apropiada: %.4f%%\n", density * 100);
        printf("   - H₁: %d elementos, H₂: %d elementos\n", h1_elements, h2_elements);
        printf("   - Lista para usar en codificación/decodificación\n");
    } else {
        printf("❌ PROBLEMAS DETECTADOS EN LA MATRIZ H\n");
        if(code.k + code.m != code.n) printf("   - Dimensiones inconsistentes\n");
        if(invalid_indices > 0) printf("   - %d índices fuera de rango\n", invalid_indices);
        if(density >= 0.01) printf("   - Densidad demasiado alta (%.4f%%)\n", density * 100);
    }

    // ⭐ DEBUG: Verificar matriz H antes de conversión
printf("\n=== DEBUG MATRIZ H ORIGINAL (COO) ===\n");

// 1) Contar grados por fila en formato COO
int *row_counts = calloc(code.m, sizeof(int));
for (int i = 0; i < code.H.num_elements; i++) {
    int row = code.H.row_indices[i];
    if (row >= 0 && row < code.m) {
        row_counts[row]++;
    } else {
        printf("❌ ERROR: row_indices[%d] = %d fuera de rango [0, %d)\n", 
               i, row, code.m);
    }
}

// 2) Mostrar estadísticas de grados
int min_deg = code.n, max_deg = 0;
float avg_deg = 0;
for (int row = 0; row < code.m; row++) {
    if (row_counts[row] < min_deg) min_deg = row_counts[row];
    if (row_counts[row] > max_deg) max_deg = row_counts[row];
    avg_deg += row_counts[row];
}
avg_deg /= code.m;

printf("Grados check nodes COO: min=%d, max=%d, avg=%.2f\n", min_deg, max_deg, avg_deg);

// 3) Mostrar filas problemáticas
printf("Filas con grado > 50:\n");
int problematic_rows = 0;
for (int row = 0; row < code.m; row++) {
    if (row_counts[row] > 50) {
        printf("  Fila %d: grado %d\n", row, row_counts[row]);
        problematic_rows++;
        if (problematic_rows >= 10) {
            printf("  ... y %d filas más\n", problematic_rows - 10);
            break;
        }
    }
}

// 4) Verificar duplicados
printf("Verificando duplicados...\n");
int duplicates = 0;
for (int i = 0; i < code.H.num_elements - 1; i++) {
    for (int j = i + 1; j < code.H.num_elements; j++) {
        if (code.H.row_indices[i] == code.H.row_indices[j] && 
            code.H.col_indices[i] == code.H.col_indices[j]) {
            duplicates++;
            if (duplicates <= 5) {
                printf("  Duplicado: (%d,%d) en índices %d y %d\n", 
                       code.H.row_indices[i], code.H.col_indices[i], i, j);
            }
        }
    }
}
printf("Total duplicados encontrados: %d\n", duplicates);

free(row_counts);
printf("=== FIN DEBUG ===\n");
    
    printf("3.1) Convirtiendo matriz H de COO a CSR...\n");

    // Extraer coordenadas COO de la matriz H
    int *coo_rows = malloc(code.H.num_elements * sizeof(int));
    int *coo_cols = malloc(code.H.num_elements * sizeof(int));

    for (int i = 0; i < code.H.num_elements; i++) {
        coo_rows[i] = code.H.row_indices[i];  // En COO: número de fila
        coo_cols[i] = code.H.col_indices[i];  // En COO: número de columna
    }

    // Convertir a CSR usando tu función existente
    SparseMatrix H_csr = coo_to_csr(code.m, code.n, code.H.num_elements, coo_rows, coo_cols);

    printf("     ✅ Conversión CSR completada: %dx%d, elementos=%d\n", 
        H_csr.num_rows, H_csr.num_cols, H_csr.num_elements);

    // Crear LDPCCode con matriz CSR
    LDPCCode code_csr = {
        .H = H_csr,
        .k = code.k,
        .n = code.n,
        .m = code.m
    };

        // ================================================================
    // ⭐ SIMULACIÓN REALISTA CON ERRORES INTRODUCIDOS
    // ================================================================

    printf("\n=== PRUEBA DECODIFICADOR MIN-SUM ===\n");

    // 1) Generar palabra código all-zeros
    int *codeword = calloc(code.n, sizeof(int));
    printf("3.1) Palabra código: all-zeros\n");

    // 2) ⭐ SNR
    float snr_db = 3.0f;  // Más alto para asegurar convergencia
    float sigma = sqrt(1.0f / (2.0f * pow(10.0f, snr_db/10.0f)));

    printf("3.2) Canal AWGN: SNR=%.1fdB, sigma=%.3f\n", snr_db, sigma);

    // 3) ⭐ SIMULACIÓN CON ERRORES FORZADOS
    float *llr_channel = malloc(code.n * sizeof(float));
    int errors_introduced = 0;

    for (int i = 0; i < code.n; i++) {
        // BPSK: bit 0 → +1, bit 1 → -1
        float transmitted = (codeword[i] == 0) ? 1.0f : -1.0f;
        
        // Ruido gaussiano fuerte
        float noise = ((float)rand()/RAND_MAX - 0.5f) * 2.0f * sigma; // ⭐ Ruido más fuerte
        float received = transmitted + noise;
        
        // ⭐ INTRODUCIR ERRORES ARTIFICIALES en ~5% de bits
        if ((rand() % 100) < 1) {
            received = -received;  // Invertir señal → forzar error
            errors_introduced++;
        }
        
        llr_channel[i] = 2.0f * received / (sigma * sigma);
    }

    printf("3.3) Errores artificiales introducidos: %d/16200 (%.2f%%)\n", 
        errors_introduced, (float)errors_introduced/code.n * 100);

    // 4) Configuración más conservadora
    DecoderConfig config = {
        .max_iter = 100,
        .alpha = 0.70f,    // Menos agresivo
        .beta = 0.30f      // Offset medio
    };

    printf("3.4) Config Min-Sum: max_iter=%d, alpha=%.2f, beta=%.2f\n", 
        config.max_iter, config.alpha, config.beta);

    // 5) ⭐ MOSTRAR LLR INICIALES para debug
    printf("3.5) DEBUG LLR iniciales [0:9]: ");
    for (int i = 0; i < 10; i++) {
        printf("%.1f ", llr_channel[i]);
    }
    printf("\n");

    // 6) Decodificar
    int *decoded_bits = malloc(code.n * sizeof(int));

    printf("3.6) 🔄 Iniciando decodificación Min-Sum...\n");
    int iterations = decode_minsum(&code, llr_channel, decoded_bits, config);

    // 7) Verificar resultado detallado
    int bit_errors = 0;
    for (int i = 0; i < code.n; i++) {
        if (decoded_bits[i] != codeword[i]) bit_errors++;
    }

    printf("3.7) 📊 Resultado: %d iteraciones, errores_finales=%d/16200 (%.3f%%)\n", 
        iterations, bit_errors, (float)bit_errors/code.n * 100);

    if (bit_errors == 0) {
        printf("     ✅ DECODIFICACIÓN PERFECTA (corrigió %d errores)\n", errors_introduced);
    } else if (bit_errors < errors_introduced) {
        printf("     ⚠️ DECODIFICACIÓN PARCIAL (corrigió %d de %d errores)\n", 
            errors_introduced - bit_errors, errors_introduced);
    } else {
        printf("     ❌ DECODIFICACIÓN FALLIDA (empeoró los errores)\n");
    }

    // 8) Debug comparativo
    printf("3.8) DEBUG comparativo:\n");
    printf("     Original [0:9]: ");
    for (int i = 0; i < 10; i++) printf("%d", codeword[i]);
    printf("\n     Decodific[0:9]: ");
    for (int i = 0; i < 10; i++) printf("%d", decoded_bits[i]);
    printf("\n");

    // 9) Limpiar
    free(codeword);
    free(llr_channel);
    free(decoded_bits);



    // 4) libera H
    freeSparseMatrix(&H);
    printf("4) fin main, H liberada\n"); fflush(stdout);
    return 0;
}
