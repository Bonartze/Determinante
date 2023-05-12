
#include <stdlib.h>
#include <stdio.h>

void fill_matrix(double *mtr, int N) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            mtr[i * N + j] = rand() % 40;
}

void printf_matrix(double *mtr, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%.2f ", mtr[i * N + j]);
        printf("\n");
    }
    printf("\n");
}