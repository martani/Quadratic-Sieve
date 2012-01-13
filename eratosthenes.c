/*
 * eratosthenes.c
 *
 *  Created on: Jan 7, 2012
 *      Author: martani
 */
/* eratosthenes sieve to calculate primes under a certain Base
 */

#include <stdio.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <gmp.h>

#define GET_BIT_AT(index) ((numbers[index>>3] & (1<<(index&7))) >> (index&7))
#define SET_BIT_AT(index) (numbers[index>>3] |= (1 << (index&7)))
#define CLEAR_BIT_AT(index) (numbers[index>>3] &= ~(1 << (index&7)))

char *numbers;
int64_t base_ref;

//Can handle up to base 2^31 * 8 * 8
int64_t sieve_primes_up_to(int64_t base) {
	int64_t num_primes = 0;
	int64_t i;

	//Find primes, base included
	base++;

	numbers = (char*) calloc(base / 64 + 1, sizeof(uint64_t));
	base_ref = base;

	for (i = 0; i < base; i++)
		SET_BIT_AT(i);

	int64_t p = 2;
	num_primes++;
	int64_t offset;

	while (1) {
		offset = 2 * p;

		while (offset < base) {
			CLEAR_BIT_AT(offset);
			offset += p;
		}

		offset = p + 1;
		while (offset < base && (GET_BIT_AT(offset) == 0)) {
			offset++;
		}

		if (offset == base)
			break;

		p = offset;
		num_primes++;
	}

	return num_primes;
}

/* Array must be malloc'ed already
 * fills the array with the first @num_primes primes calculated with sieve_primes
 */
void fill_primes(int64_t *primes_array) {
	int64_t j, i;

	for (j = 0, i = 2; i < base_ref; i++) {
		if (GET_BIT_AT(i) == 1) {
			primes_array[j] = i;
			j++;
		}
	}

	free(numbers);
}

/* Fill the array with only primes where n is a quadratic residue: x² = n (mod p) */
int fill_primes_with_quadratic_residue(int64_t *primes_array, mpz_t n) {
	int64_t j, i;
	mpz_t b;
	mpz_init(b);

	primes_array[0] = 2;
	for (j = 1, i = 3; i < base_ref; i++) {
		mpz_set_ui(b, (unsigned long) i);
		if ((GET_BIT_AT(i)) == 1 && mpz_jacobi(n, b) == 1) {
			primes_array[j] = i;
			j++;
		}
	}

	free(numbers);
	return j;
}
