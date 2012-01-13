/*
 * eratosthenes.h
 *
 *  Created on: Jan 7, 2012
 *      Author: martani
 */

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

int64_t sieve_primes_up_to(int64_t base);

/* Array must be malloc'ed already
 * fills the array with the first @num_primes primes calculated with sieve_primes
 */
void fill_primes(int64_t *primes_array);

/* Fill the array with only primes where n is a quadratic residue: x² = n (mod p) */
int fill_primes_with_quadratic_residue(int64_t *primes_array, mpz_t n);
