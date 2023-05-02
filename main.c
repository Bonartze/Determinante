#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include <string.h>

#define N 3
#define M 4

int max(int a, int b) {
    return a > b ? a : b;
}

void divide_lead_elem_0_proc(int *matrix, int rank, int total_processes, int s_c, int lead_el) {
    for (int i = 0; i < M; i++)
        matrix[s_c + i] /= matrix[s_c + lead_el];
}

void divide_leading_element_line(int *matrix, int f_l, int l_l, int l_c,
                                 int *lead_string) {
    for (int i = 0; i < f_l - l_l + 1; i++)
        for (int j = l_c; j < M; j++)
            matrix[i * N + j] -= lead_string[j] * matrix[i * N + l_c];
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
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int matrix[M * N];
    int lead_string[M];
    int lead_str = 0, lead_col = 0;
    int first_line = (N * (rank - 1)) / size;
    first_line = first_line == 0 ? first_line + 1 : first_line;
    int last_line = (N * (rank)) / size;
    int sub_matrix[last_line - first_line + 1];
    if (rank == 0) {
        fill_matrix(matrix, N, M);
        //print_matrix(matrix, N, M, rank, size);
        //   printf("\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                if (matrix[i * N + j] != 0) {
                    divide_lead_elem_0_proc(matrix, rank, size, i, j);
                    memcpy(lead_string, matrix, sizeof(int) * M);
                    for (int process = 1; process < size; process++) {
                        first_line = (N * (process - 1)) / size;
                        last_line = (N * (process)) / size;
                        int ind = 0;
                        for (int k = first_line; k < last_line; k++)
                            for (int t = 0; t < M; t++)
                                sub_matrix[ind++] = matrix[i * N + j];
                        MPI_Send(sub_matrix, last_line - first_line + 1, MPI_INT, process, MPI_ANY_TAG, MPI_COMM_WORLD);
                        MPI_Send(lead_string, M, MPI_INT, process, MPI_ANY_TAG, MPI_COMM_WORLD);
                        MPI_Send(&i, 1, MPI_INT, process, MPI_ANY_TAG, MPI_COMM_WORLD);
                        MPI_Recv(sub_matrix, last_line - first_line + 1, MPI_INT, process, MPI_ANY_TAG, MPI_COMM_WORLD,
                                 0);
                        int temp = 0;
                        for (int a = first_line; a < last_line; a++)
                            for (int b = 0; b < M; b++)
                                matrix[a * N + b] = sub_matrix[temp++];
                    }
                    break;
                }
            }
        }
        //  print_matrix(matrix, N, M, rank, size);
    } else {
        for (int i = 0; i < N; i++) {
            int num_lead_str;
            MPI_Recv(sub_matrix, last_line - first_line + 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
            MPI_Recv(lead_string, M, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
            MPI_Recv(&num_lead_str, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
            divide_leading_element_line(sub_matrix, first_line, last_line, num_lead_str, lead_string);
            MPI_Send(sub_matrix, last_line - first_line + 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
}
