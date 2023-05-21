#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include <string.h>
#include "matrix.h" //for using operations with the matrices

#define N 4

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n_rows = N / size;   //number of rows in one process
    const int start_row = n_rows * rank;
    const int finish_row = start_row + n_rows;
    double mtr[N * N];
    double m_chunk[N * n_rows];
    double pivot_row[N];
    if (rank == 0) {
        fill_matrix(mtr, N);
        printf_matrix(mtr, N);
    }
    MPI_Scatter(mtr, N * n_rows, MPI_DOUBLE, m_chunk, N * n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Request requests[size];
    for (int row = 0; row < finish_row; row++) {
        int mapped_rank = row / n_rows;
        if (rank == mapped_rank) {
            int local_row = row % n_rows;
            double pivot = m_chunk[local_row * N + row];
            for (int i = mapped_rank + 1; i < size; i++) {
                MPI_Isend(m_chunk + N * local_row, N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
                          &requests[i]);  //sending current row
            }
            for (int elim_row = local_row + 1; elim_row < n_rows; elim_row++) {
                double scale = m_chunk[elim_row * N + row];
                for (int col = row; col < N; col++) {
                    m_chunk[elim_row * N + col] -= ((m_chunk[local_row * N + col] * scale) / pivot);
                }
            }
            for (int i = mapped_rank + 1; i < size; i++) {
                MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
            }
        } else {
            MPI_Recv(pivot_row, N, MPI_DOUBLE, mapped_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int elim_row = 0; elim_row < n_rows; elim_row++) {
                double scale = m_chunk[elim_row * N + row];
                for (int col = row; col < N; col++) {
                    m_chunk[elim_row * N + col] -= pivot_row[col] * scale / pivot_row[row];
                }
            }
        }
    }
    MPI_Gather(m_chunk, n_rows * N, MPI_DOUBLE, mtr, n_rows * N,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf_matrix(mtr, N);
    }
    MPI_Finalize();
    return 0;
}