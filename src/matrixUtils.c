#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "matrixUtils.h"

int* addcr(const int *c, const int *r, int M, int N) {
    int *A = (int *)malloc(M * N * sizeof(int));
    if (A == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            A[m * N + n] = r[n] + c[m];
        }
    }

    return A;
}

SparseMatrixCSR* coo_to_csr(SparseMatrix *coo){
    int nrows = coo->num_rows;
    int nnz = coo->num_elements;

    SparseMatrixCSR *csr = malloc(sizeof(SparseMatrixCSR));
    csr->num_rows = nrows;
    csr->num_cols = coo->num_cols;
    csr->num_elements = nnz;
    csr->indptr = calloc(nrows + 1, sizeof(int));
    csr->indices = malloc(nnz * sizeof(int));
    csr->values = malloc(nnz * sizeof(int));

    if (!csr->indptr || !csr->indices || !csr->values) {
        free(csr); return NULL;
    }

    // Paso 1: contar los elementos por fila
    for(int i = 0; i < nnz; i++) {
        csr->indptr[coo->row_indices[i] + 1]++;
    }

    // Paso 2: acumulado para indptr
    for(int i = 1; i <= nrows; i++) {
        csr->indptr[i] += csr->indptr[i - 1];
    }

    // Copiar pos inicial
    int *temp = malloc((nrows + 1) * sizeof(int));
    memcpy(temp, csr->indptr, (nrows + 1) * sizeof(int));

    // Paso 3: llenar índices y valores
    for(int i = 0; i < nnz; i++) {
        int r = coo->row_indices[i];
        int dest = temp[r]++;
        csr->indices[dest] = coo->col_indices[i];
        csr->values[dest] = coo->values[i];
    }

    free(temp);
    return csr;
}




int* matrixToVector(int **matrix, int rows, int cols){

	int *vec = (int*)malloc(rows * cols * sizeof(int));
    if (vec == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    int idx = 0;
    for (int col = 0; col < cols; col++) {
        for (int row = 0; row < rows; row++) {
            vec[idx++] = matrix[col][row];
        }
    }

    return vec;
}

void freeMatrix(int **matrix, int rows) {
    if (matrix == NULL) return;
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void printMatrix(int **matrix, int rows, int cols){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int isSystematicForm(SparseMatrix H, int k) {
    int m = H.num_rows;
    int n = H.num_cols;

    if (n - k != m) {
        printf("ERROR: Para forma sistemática, m debe ser igual a n - k\n");
        return 0;
    }

    for (int row = 0; row < m; row++) {
        int found_one = 0;
        for (int i = H.row_indices[row]; i < H.row_indices[row + 1]; i++) {
            int col = H.col_indices[i];
            if (col < k) continue;

            if (col - k == row && !found_one) {
                found_one = 1;
            } else {
                return 0;
            }
        }
        if (!found_one) return 0;
    }

    return 1;
}

/*SparseLIL csr_to_lil(const SparseMatrix *csr) {
    SparseLIL lil;
    lil.num_rows = csr->num_rows;
    lil.num_cols = csr->num_cols;

    lil.cols = malloc(lil.num_rows * sizeof(int*));
    lil.counts = malloc(lil.num_rows * sizeof(int));

    if (!lil.cols || !lil.counts) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    for (int row = 0; row < lil.num_rows; row++) {
        int start = csr->row_indices[row];
        int end   = csr->row_indices[row + 1];
        int count = end - start;

        lil.counts[row] = count;
        lil.cols[row] = malloc(count * sizeof(int));
        if (!lil.cols[row]) {
            perror("malloc fila lil");
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < count; i++) {
            lil.cols[row][i] = csr->col_indices[start + i];
        }
    }

    return lil;
}

void free_lil(SparseLIL *lil) {
    for (int i = 0; i < lil->num_rows; i++) {
        free(lil->cols[i]);
    }
    free(lil->cols);
    free(lil->counts);
    lil->cols = NULL;
    lil->counts = NULL;
}

int gaussian_elimination_systematic(SparseLIL *lil, int k, int *col_perm) {
    int m = lil->num_rows;
    int n = lil->num_cols;

    if (n - k != m) {
        printf("ERROR: n - k != m, no se puede sistematizar\n");
        return 0;
    }

    // Inicializar permutación como identidad
    for (int i = 0; i < n; i++) {
        col_perm[i] = i;
    }

    for (int row = 0; row < m; row++) {
        int target_col = k + row;
        int pivot_index = -1;

        // Buscar si ya existe el 1 en la columna deseada
        for (int i = 0; i < lil->counts[row]; i++) {
            if (lil->cols[row][i] == target_col) {
                pivot_index = i;
                break;
            }
        }

        // Si no está, intentar encontrar otro 1 en la fila y permutar
        if (pivot_index == -1) {
            int swap_col = -1;
            for (int i = 0; i < lil->counts[row]; i++) {
                int col = lil->cols[row][i];
                if (col >= k && col < n) {
                    swap_col = col;
                    break;
                }
            }

            if (swap_col == -1) {
                printf("❌ No se puede sistematizar, fila %d no tiene pivot válido\n", row);
                return 0;
            }

            // Permutamos columnas swap_col <-> target_col
            for (int r = 0; r < m; r++) {
                for (int i = 0; i < lil->counts[r]; i++) {
                    if (lil->cols[r][i] == swap_col) {
                        lil->cols[r][i] = target_col;
                    } else if (lil->cols[r][i] == target_col) {
                        lil->cols[r][i] = swap_col;
                    }
                }
            }

            // Actualizar permutación
            int tmp = col_perm[swap_col];
            col_perm[swap_col] = col_perm[target_col];
            col_perm[target_col] = tmp;
        }

        // Hacer cero el resto de esa columna (XOR entre filas)
        for (int r = 0; r < m; r++) {
            if (r == row) continue;

            // ¿Tiene un 1 en target_col?
            int has_one = 0;
            for (int i = 0; i < lil->counts[r]; i++) {
                if (lil->cols[r][i] == target_col) {
                    has_one = 1;
                    break;
                }
            }

            if (has_one) {
                xor_rows(lil, r, row);  // fila_r = fila_r XOR fila_row
            }
        }
    }

    printf("✅ Matriz H reducida a forma sistemática mediante XOR + permutación\n");
    return 1;
}

void xor_rows(SparseLIL *lil, int dest_row, int src_row) {
    int *a = lil->cols[dest_row];
    int *b = lil->cols[src_row];
    int na = lil->counts[dest_row];
    int nb = lil->counts[src_row];

    // Resultado temporal (union simétrica)
    int *temp = malloc((na + nb) * sizeof(int));
    int nt = 0;

    int i = 0, j = 0;
    while (i < na && j < nb) {
        if (a[i] < b[j]) {
            temp[nt++] = a[i++];
        } else if (a[i] > b[j]) {
            temp[nt++] = b[j++];
        } else {
            // Mismo índice en ambas → 1 XOR 1 = 0, se elimina
            i++;
            j++;
        }
    }

    while (i < na) temp[nt++] = a[i++];
    while (j < nb) temp[nt++] = b[j++];

    free(lil->cols[dest_row]);
    lil->cols[dest_row] = temp;
    lil->counts[dest_row] = nt;
}

SparseMatrix lil_to_csr(const SparseLIL *lil) {
    SparseMatrix csr;
    csr.num_rows = lil->num_rows;
    csr.num_cols = lil->num_cols;

    // Calcular número total de elementos no nulos
    int nnz = 0;
    for (int i = 0; i < lil->num_rows; i++) {
        nnz += lil->counts[i];
    }
    csr.num_elements = nnz;

    csr.row_indices = calloc(lil->num_rows + 1, sizeof(int));
    csr.col_indices = malloc(nnz * sizeof(int));
    csr.values = NULL; // los valores en GF(2) se asumen 1

    if (!csr.row_indices || !csr.col_indices) {
        perror("malloc CSR");
        exit(EXIT_FAILURE);
    }

    int index = 0;
    for (int i = 0; i < lil->num_rows; i++) {
        csr.row_indices[i] = index;
        for (int j = 0; j < lil->counts[i]; j++) {
            csr.col_indices[index++] = lil->cols[i][j];
        }
    }
    csr.row_indices[lil->num_rows] = index;

    return csr;
}*/

uint8_t **csr_to_dense(const SparseMatrix *csr) {
    int m = csr->num_rows;
    int n = csr->num_cols;

    // Reservamos matriz densa de m × n e inicializamos a 0
    uint8_t **dense = malloc(m * sizeof(uint8_t*));
    for (int i = 0; i < m; i++) {
        dense[i] = calloc(n, sizeof(uint8_t));  // calloc para rellenar con ceros
    }

    // Rellenamos con 1 donde haya valores en la CSR
    for (int row = 0; row < m; row++) {
        int start = csr->row_indices[row];
        int end   = csr->row_indices[row + 1];
        for (int i = start; i < end; i++) {
            int col = csr->col_indices[i];
            dense[row][col] = 1;  // valores siempre son 1 en GF(2)
        }
    }

    return dense;
}

void free_dense(uint8_t **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

uint8_t **compute_nullspace_gf2(uint8_t **H, int m, int n, int *k_out) {
    int rows = m;
    int cols = n;

    // Matriz aumentada [H | I]
    int aug_cols = cols + rows;
    uint8_t **A = malloc(rows * sizeof(uint8_t*));
    for (int i = 0; i < rows; i++) {
        A[i] = calloc(aug_cols, sizeof(uint8_t));
        for (int j = 0; j < cols; j++) {
            A[i][j] = H[i][j];
        }
        A[i][cols + i] = 1;  // matriz identidad añadida a la derecha
    }

    // Eliminación gaussiana en GF(2)
    int *pivot_col = malloc(rows * sizeof(int));
    for (int i = 0; i < rows; i++) pivot_col[i] = -1;

    int pivot_row = 0;
    for (int col = 0; col < cols && pivot_row < rows; col++) {
        int sel = -1;
        for (int r = pivot_row; r < rows; r++) {
            if (A[r][col]) {
                sel = r;
                break;
            }
        }

        if (sel == -1) continue;

        // intercambiar filas si hace falta
        if (sel != pivot_row) {
            uint8_t *tmp = A[sel];
            A[sel] = A[pivot_row];
            A[pivot_row] = tmp;
        }

        pivot_col[pivot_row] = col;

        // hacer ceros arriba y abajo
        for (int r = 0; r < rows; r++) {
            if (r != pivot_row && A[r][col]) {
                for (int j = 0; j < aug_cols; j++) {
                    A[r][j] ^= A[pivot_row][j];
                }
            }
        }

        pivot_row++;
    }

    int piv = pivot_row;
    int k = cols - piv;
    *k_out = k;

    // Almacenar índices de columnas libres
    int *free_cols = malloc(k * sizeof(int));
    int fc = 0;
    int current_pivot = 0;

    for (int col = 0; col < cols; col++) {
        if (current_pivot < piv && pivot_col[current_pivot] == col) {
            current_pivot++;
        } else {
            free_cols[fc++] = col;
        }
    }

    // Construir G: una fila por columna libre
    uint8_t **G = malloc(k * sizeof(uint8_t*));
    for (int i = 0; i < k; i++) {
        G[i] = calloc(cols, sizeof(uint8_t));
        G[i][free_cols[i]] = 1;

        // recorrer filas de la matriz reducida
        for (int r = 0; r < piv; r++) {
            int pc = pivot_col[r];
            if (A[r][free_cols[i]]) {
                G[i][pc] = 1;
            }
        }
    }

    // Liberar memoria
    for (int i = 0; i < rows; i++) free(A[i]);
    free(A);
    free(free_cols);
    free(pivot_col);

    return G;
}


// Función opcional para validar H·Gᵗ = 0
int check_HGt(uint8_t **H, int m, int n, uint8_t **G, int k) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            int sum = 0;
            for (int l = 0; l < n; l++) {
                sum ^= H[i][l] & G[j][l];
            }
            if (sum != 0) {
                printf("❌ H·Gᵗ != 0 en fila %d columna %d\n", i, j);
                return 0;
            }
        }
    }
    return 1;
}

int gcd(int a, int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        int r = a % b;
        a = b;
        b = r;
    }
    return a;
}

int saveG(const char *filename, uint8_t **G, int k, int n){
    FILE *f = fopen(filename, "wb");
    if (!f) return -1;
    if (fwrite(&k, sizeof(int), 1, f)!=1 ||
        fwrite(&n, sizeof(int), 1, f)!=1) {
        fclose(f);
        return -1;
    }
    for (int i = 0; i < k; i++) {
        if (fwrite(G[i], sizeof(uint8_t), n, f) != (size_t)n) {
            fclose(f);
            return -1;
        }
    }
    fclose(f);
    return 0;
}

uint8_t **loadG(const char *filename, int *k_out, int *n_out){
    struct stat st;
    if (stat(filename, &st)!=0) return NULL;
    FILE *f = fopen(filename, "rb");
    if (!f) return NULL;
    int k,n;
    if (fread(&k, sizeof(int), 1, f)!=1 ||
        fread(&n, sizeof(int), 1, f)!=1) {
        fclose(f);
        return NULL;
    }
    uint8_t **G = malloc(k * sizeof(uint8_t*));
    if (!G) { fclose(f); return NULL; }
    for (int i = 0; i < k; i++) {
        G[i] = malloc(n * sizeof(uint8_t));
        if (fread(G[i], sizeof(uint8_t), n, f) != (size_t)n) {
            for (int j = 0; j <= i; j++) free(G[j]);
            free(G);
            fclose(f);
            return NULL;
        }
    }
    fclose(f);
    *k_out = k;
    *n_out = n;
    return G;
}