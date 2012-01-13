/*
 * custom_types.h
 *
 *  Created on: Jan 12, 2012
 *      Author: martani
 */

#include <gmp.h>

#ifndef CUSTOM_TYPES_H_
#define CUSTOM_TYPES_H_

typedef struct {
	mpz_t value_x;
	mpz_t value_x_squared;
	/* uint64_t *factors_exp;     /* this can be used to keep track for the full factorization of
                                   * an element x²-n. For the use of this version, refer to the
                                   * qs_exponent_vector.c version
                                   */

	mpz_t factors_vect; /* this nb_primes_in_base bit vector is a substitute of the
                         * factors_exp array, all we actually need for the elimination
                         * in the Linear algebra step is the exponents modulo 2.
                         * This saves huge space within the sieving step
                         */
} smooth_number_t;

typedef struct modular_root{
  unsigned long root1;
  unsigned long root2;
} modular_root_t;



#endif
