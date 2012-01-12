/*
 * matrix.c
 *
 *  Created on: Jan 7, 2012
 *      Author: martani
 */

#include <gmp.h>
#include <inttypes.h>

#ifndef MATRIX_H_
#define MATRIX_H_

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

typedef struct {
	mpz_t *MATRIX;
	mpz_t *IDENTITY;
	uint64_t rows;
	uint64_t cols;
	uint64_t next_free_row;		/* used to insert rows in the matrix, point to the next free row to be inserted position*/
} matrix_t;

/* allocates space for m rows * n columns matrix for MATRIX and IDENTITY */
void init_matrix(matrix_t *matrix, uint64_t m, uint64_t n);

void push_row(matrix_t *matrix, mpz_t row);

void print_matrix_matrix(matrix_t *matrix);

void print_matrix_identity(matrix_t *matrix);

/* performs a Gauss elimination on matrix->MATRIX, result (linear dependence) will be in the matrix->IDENTITY */
void gauss_elimination(matrix_t *matrix);

/* does not check for bounds, the caller must */
void get_matrix_row(mpz_t rop, matrix_t *matrix, uint64_t row_index);

void get_identity_row(mpz_t rop, matrix_t *matrix, uint64_t row_index);

int test();

#endif /* MATRIX_H_ */
