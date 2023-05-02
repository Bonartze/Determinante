#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>

#define N 3
#define M 4

int max(int a, int b) {
    return a > b ? a : b;
}

void divide_lead_elem_0_proc(int *matrix, int rank, int total_processes, int s_c, int lead_el) {
    for (int i = 0; i < M; i++)
        matrix[s_c + i] /= matrix[s_c + lead_el];
}

void divide_leading_element_line(int *matrix, int rank, int total_processes, int f_l, int l_l, int l_s, int l_c) {
    for (int i = max(f_l, l_s); i < l_l; i++)
        for (int j = l_c; j < M;
             j++)
            matrix[i * N + j] -= matrix[l_s * N + l_c] * matrix[i * N + l_c];

}

void fill_matrix(int *mtr, int n, int m) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            mtr[i * n + j] = rand() % 30;
}

void print_matrix(int *mtr, int n, int m, int rank, int total_proc) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            printf("%i ", mtr[i * n + j]);
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    int matrix[M * N];
    int size, rank;
    int lead_str = 0, lead_col = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int first_line = (N * (rank - 1)) / size;
    first_line = first_line == 0 ? first_line + 1 : first_line;
    int last_line = (N * (rank)) / size;
    if (rank == 0) {
        fill_matrix(matrix, N, M);
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < M; i++) {
                if (matrix[j * N + i] != 0) {
                    divide_lead_elem_0_proc(matrix, rank, size, lead_str, lead_col);
                    MPI_Bcast(matrix, N * M, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(&j, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    break;
                }
            }
        }
        print_matrix(matrix, N, M, rank, size);
    } else {
        for (int i = 1; i < N; i++) {
            printf("%i\n", i);
            MPI_Bcast(matrix, N * M, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&lead_col, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&lead_str, 1, MPI_INT, 0, MPI_COMM_WORLD);
            divide_leading_element_line(matrix, rank, size, first_line, lead_str, last_line, lead_col);
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
}
