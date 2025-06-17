#include "dvbs2ldpcShort.h"
#include "matrixUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

    // 4) libera H
    freeSparseMatrix(&H);
    printf("4) fin main, H liberada\n"); fflush(stdout);
    return 0;
}
