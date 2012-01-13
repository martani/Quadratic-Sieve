/*
 * matrix.c
 *
 *  Created on: Jan 7, 2012
 *      Author: martani
 */
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <inttypes.h>

#include "matrix.h"

/* allocates space for m rows * n columns matrix for MATRIX and IDENTITY */
void init_matrix(matrix_t *matrix, uint64_t m, uint64_t n) {
	matrix->MATRIX = (mpz_t *) calloc(m, sizeof(mpz_t));
	matrix->IDENTITY = (mpz_t *) calloc(m, sizeof(mpz_t));
	matrix->rows = m;
	matrix->cols = n;
	matrix->next_free_row = 0;
}

void push_row(matrix_t *matrix, mpz_t row) {
	mpz_set(matrix->MATRIX[matrix->next_free_row], row);

	mpz_init2(matrix->IDENTITY[matrix->next_free_row], matrix->cols); /* initializes a n bit vector all set to 0 */
	mpz_setbit(matrix->IDENTITY[matrix->next_free_row], matrix->next_free_row); /* set the next_free_row bit to 1 */
	matrix->next_free_row++;
}

void print_matrix_matrix(matrix_t *matrix) {
	uint64_t i, j;

	printf("\nMATRIX\n");
	for (i = 0; i < matrix->rows; i++) {
		printf("[");
		for (j = 0; j < matrix->cols; j++) {
			printf("%2d", mpz_tstbit(matrix->MATRIX[i], j));
		}
		printf(" ]\n");
	}
	printf("\n");
}

void print_matrix_identity(matrix_t *matrix) {
	uint64_t i, j;

	printf("\nIDENTITY\n");
	for (i = 0; i < matrix->rows; i++) {
		printf("[");
		for (j = 0; j < matrix->cols; j++) {
			printf("%2d", mpz_tstbit(matrix->IDENTITY[i], j));
		}
		printf(" ]\n");
	}
	printf("\n");
}

void free_matrix(matrix_t *matrix) {
	free(matrix->MATRIX);
	free(matrix->IDENTITY);
}

/* performs a Gauss elimination on matrix->MATRIX, result (linear dependence) will be in the matrix->IDENTITY */
void gauss_elimination(matrix_t *matrix) {
	printf("\nPerforming Gauss elimination..\n");
	mpz_t *m = matrix->MATRIX;
	mpz_t *I = matrix->IDENTITY;

	uint64_t col, row, next_row, next_pivot;
	for (next_row = 0, col = 0; col < min(matrix->cols, matrix->rows); col++) /* for all rows*/
	{
		next_pivot = -1;
		for (row = next_row; row < matrix->rows; row++) /* search for the next pivot*/
		{
			if (mpz_tstbit(m[row], col)) {
				next_pivot = row; /* row contains the next pivot */
				next_row++;
				break;
			}
		}

		if (next_pivot == -1)
			continue;

		if (next_pivot != next_row - 1) /* current row is not the pivot, switch rows */
		{
			mpz_swap(m[next_pivot], m[next_row - 1]);
			mpz_swap(I[next_pivot], I[next_row - 1]);
		}

		for (row = next_row; row < matrix->rows; row++) {
			if (mpz_tstbit(m[row], col)) {
				mpz_xor(m[row], m[row], m[next_row - 1]); /* XOR the rows to eliminate the 1 in position (row, next_row-1)*/
				mpz_xor(I[row], I[row], I[next_row - 1]);
			}
		}
	}
}

/* does not check for bounds, the caller must */
void get_matrix_row(mpz_t rop, matrix_t *matrix, uint64_t row_index) {
	mpz_set(rop, matrix->MATRIX[row_index]);
}

void get_identity_row(mpz_t rop, matrix_t *matrix, uint64_t row_index) {
	mpz_set(rop, matrix->IDENTITY[row_index]);
}

int test() {
	matrix_t matrix;

	init_matrix(&matrix, 6, 6);
	mpz_t r1, r2, r3;
	mpz_init(r1);
	mpz_init(r2);
	mpz_init(r3);

	mpz_set_ui(r1, 12);
	mpz_set_ui(r2, 11);
	mpz_set_ui(r3, 23);

	push_row(&matrix, r1);
	push_row(&matrix, r2);
	push_row(&matrix, r3);

	mpz_set_ui(r3, 1);
	push_row(&matrix, r3);

	mpz_set_ui(r3, 1);
	push_row(&matrix, r3);

	mpz_set_ui(r3, 23);
	push_row(&matrix, r3);

	print_matrix_matrix(&matrix);
	print_matrix_identity(&matrix);

	gauss_elimination(&matrix);

	print_matrix_matrix(&matrix);
	print_matrix_identity(&matrix);

	return 0;
}
