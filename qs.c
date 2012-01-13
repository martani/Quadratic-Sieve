/*
 * qs.c
 *
 *  Created on: Dec 25, 2011
 *      Author: martani
 */
/*-------------- Factorize numbers with the Quadratic Sieve method --------------*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <inttypes.h>
#include <sys/resource.h>
#include <sys/time.h>

#include "eratosthenes.h"
#include "gmp_patch.h"
#include "matrix.h"

struct timeval start_global;
struct timeval end_global;
struct timeval start;
struct timeval end;
struct timeval elapsed;

typedef struct modular_root {
	unsigned long root1;
	unsigned long root2;
} modular_root_t;

typedef struct {
	mpz_t value_x;
	mpz_t value_x_squared;
	/* uint64_t *factors_exp;      /* this can be used to keep track for the full factorization of
	 * an element x²-n. For the use of this version, refer to the
	 * qs_exponent_vector.c version
	 */

	mpz_t factors_vect; /* this nb_primes_in_base bit vector is a substitute of the
	 * factors_exp array, all we actually need for the elimination
	 * in the Linear algebra step is the exponents modulo 2.
	 * This saves huge space within the sieving step
	 */
} smooth_number_t;

mpz_t N; /* number to factorize */
matrix_t matrix;

uint64_t nb_smooth_numbers_found = 0;
uint64_t nb_qr_primes; /* number of primes p where N is a quadratic residue mod p */
uint64_t *primes; /* array holding the primes of the smoothness base */

smooth_number_t *smooth_numbers;
int NB_VECTORS_OFFSET = 5; /* number of additional rows in the matrix, to make sure that a linear relation exists */

//-----------------------------------------------------------
// base <- exp((1/2) sqrt(ln(n) ln(ln(n))))
//-----------------------------------------------------------
void get_smoothness_base(mpz_t base, mpz_t n) {
	mpfr_t fN, lnN, lnlnN;
	mpfr_init(fN), mpfr_init(lnN), mpfr_init(lnlnN);

	mpfr_set_z(fN, n, MPFR_RNDU);
	mpfr_log(lnN, fN, MPFR_RNDU);
	mpfr_log(lnlnN, lnN, MPFR_RNDU);

	mpfr_mul(fN, lnN, lnlnN, MPFR_RNDU);
	mpfr_sqrt(fN, fN, MPFR_RNDU);
	mpfr_div_ui(fN, fN, 2, MPFR_RNDU);
	mpfr_exp(fN, fN, MPFR_RNDU);

	mpfr_get_z(base, fN, MPFR_RNDU);

	mpfr_clears(fN, lnN, lnlnN, NULL);
}

// returns the index of the first element to start sieving from
// (first multiple of root that is directly greater than start)
// res = p*t + root >= start */
void get_sieving_start_index(mpz_t res, mpz_t start, mpz_t p,
		unsigned long root) {
	mpz_t q, r;
	mpz_init(q);
	mpz_init(r);

	mpz_sub_ui(start, start, root);
	mpz_fdiv_qr(q, r, start, p);

	if (mpz_cmp_ui(r, 0) != 0)
		mpz_add_ui(q, q, 1);

	mpz_mul(q, q, p); /* next element p*q+root that is directly >= start */
	mpz_add_ui(q, q, root);
	mpz_set(res, q);
	mpz_clear(q);
	mpz_clear(r);
}

// given an array of exponents of factors on the base [primes], reconstructs the number
void reconstruct_mpz(mpz_t rop, uint64_t *factors_exp) {
	uint64_t i;
	mpz_t t;
	mpz_init(t);
	mpz_t p_exp;
	mpz_init(p_exp);

	mpz_set_ui(t, 1);
	for (i = 0; i < nb_qr_primes; i++) {
		mpz_set_ui(p_exp, primes[i]);
		mpz_pow_ui(p_exp, p_exp, factors_exp[i]);
		mpz_mul(t, t, p_exp);
	}

	mpz_set(rop, t);
}

// save the sooth number n to the smooth_numbers array, and at the same time its exponents vector to the matrix
mpz_t tmp_matrix_row;
void save_smooth_number(smooth_number_t n) {
	mpz_clear(tmp_matrix_row); /* tmp_matrix_row must be initialized already */
	mpz_init2(tmp_matrix_row, nb_qr_primes); /* init a vector of *exactly* nb_qr_primes bits */

	if (nb_smooth_numbers_found > nb_qr_primes + NB_VECTORS_OFFSET - 1) /* if we have sufficient smooth numbers, skip saving */
		return;

	smooth_number_t tmp;
	mpz_init(tmp.value_x);
	mpz_init(tmp.value_x_squared);

	mpz_set(tmp.value_x, n.value_x);
	mpz_pow_ui(tmp.value_x_squared, n.value_x, 2);
	mpz_sub(tmp.value_x_squared, tmp.value_x_squared, N); /* x²-N. Saving this will enable us to not go through exponents
	 * and reconstruct the original number */

	/* otherwise we can reconstruct value_x_squared from the exponents vector, this is useful in the factoring step
	 * to calculate the square modulo N from the factors directly.
	 * It takes a lot of space, doesn't woth it anyways */

	//tmp.factors_exp = calloc(nb_qr_primes, sizeof(uint64_t));
	//memcpy(tmp.factors_exp, n.factors_exp, nb_qr_primes * sizeof(uint64_t));
	//reconstruct_mpz(tmp.value_x_squared, tmp.factors_exp);
	smooth_numbers[nb_smooth_numbers_found++] = tmp;

	/*** reconstruct and saves the smooth number to the GF2 matrix ***//* already done in the sieving step */
	/* uint64_t i;
	 for(i=0; i<nb_qr_primes; i++)
	 {
	 if(n.factors_exp[i]&1)
	 mpz_setbit(tmp_matrix_row, i);
	 }*/

	/* the coefficient vector in GF2 has already been constructed */
	push_row(&matrix, n.factors_vect);
}

void print_lib_version() {
	printf("GMP : %s\n", gmp_version);
	printf("MPFR library: %-12s\n\n", mpfr_get_version());
}

// seems to not work properly if user doesn't have enough privileges
void show_mem_usage() {
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	printf("\nUsing: %.2fMB of memory\n", ru.ru_maxrss / 1024.0);
}

void START_TIMER() {
	gettimeofday(&start, NULL);
}

void STOP_TIMER_PRINT_TIME(char *s) {
	gettimeofday(&end, NULL);
	timersub(&end, &start, &elapsed);
	printf("%s [Time: %.3f ms]\n", s,
			elapsed.tv_sec * 1000 + elapsed.tv_usec / (double) 1000);
}

FILE *fp;
void OPEN_LOG_FILE(char *file_name) {
	char s[512];
	sprintf(s, "%s_%s_exp", file_name, mpz_get_str(NULL, 10, N));
	if ((fp = fopen(s, "w")) == NULL) {
		printf("Cannot open file.\n");
		fp = NULL;
	}
}

void APPEND_TO_LOG_FILE(char *text) {
	if (fp != NULL)
		fprintf(fp, "%s\n", text);
}

int main(int argc, char **argv) {
	gettimeofday(&start_global, NULL);
	print_lib_version();

	mpz_init(N);
	mpz_t B;
	mpz_init(B);

	unsigned long int uBase;
	int64_t nb_primes;
	modular_root_t *modular_roots;

	uint64_t i, j;

	if (mpz_init_set_str(N, argv[1], 10) == -1) {
		printf("Cannot load N %s\n", argv[1]);
		exit(2);
	}

	mpz_t sqrtN, rem;
	mpz_init(sqrtN);
	mpz_init(rem);
	mpz_sqrtrem(sqrtN, rem, N);

	if (mpz_cmp_ui(rem, 0) != 0) /* if not perfect square, calculate the ceiling */
		mpz_add_ui(sqrtN, sqrtN, 1);
	else /* N is a perfect square, factored! */
	{
		printf("\n<<<[FACTOR]>>> %s\n", mpz_get_str(NULL, 10, sqrtN));
		return 0;
	}

	if (mpz_probab_prime_p(N, 10) > 0) /* don't bother factoring */
	{
		printf("N:%s is prime\n", mpz_get_str(NULL, 10, N));
		exit(0);
	}

	OPEN_LOG_FILE("freq");

//--------------------------------------------------------
//  calculate the smoothness base for the given N
//--------------------------------------------------------
	get_smoothness_base(B, N); /* if N is too small, the program will surely fail, please consider a pen and paper instead */
	uBase = mpz_get_ui(B);
	printf("n: %s\tBase: %s\n", mpz_get_str(NULL, 10, N),
			mpz_get_str(NULL, 10, B));

//--------------------------------------------------------
// sieve primes that are less than the smoothness base using Eratosthenes sieve
//--------------------------------------------------------
	START_TIMER();
	nb_primes = sieve_primes_up_to((int64_t) (uBase));

	printf("\nPrimes found %" PRId64 " [Smoothness Base %lu]\n", nb_primes,
			uBase);
	STOP_TIMER_PRINT_TIME("\tEratosthenes Sieving done");

//--------------------------------------------------------
// fill the primes array with primes to which n is a quadratic residue
//--------------------------------------------------------
	START_TIMER();
	primes = calloc(nb_primes, sizeof(int64_t));
	nb_qr_primes = fill_primes_with_quadratic_residue(primes, N);

	/*for(i=0; i<nb_qr_primes; i++)
	 printf("%" PRId64 "\n", primes[i]);*/

	printf("\nN-Quadratic primes found %" PRId64 "\n", nb_qr_primes);
	STOP_TIMER_PRINT_TIME("\tQuadratic prime filtering done");

//--------------------------------------------------------
// calculate modular roots
//--------------------------------------------------------
	START_TIMER();
	modular_roots = calloc(nb_qr_primes, sizeof(modular_root_t));
	mpz_t tmp, r1, r2;
	mpz_init(tmp);
	mpz_init(r1);
	mpz_init(r2);

	for (i = 0; i < nb_qr_primes; i++) {
		mpz_set_ui(tmp, (unsigned long) primes[i]);
		mpz_sqrtm(r1, N, tmp); /* calculate the modular root */
		mpz_neg(r2, r1); /* -q mod n */
		mpz_mod(r2, r2, tmp);

		modular_roots[i].root1 = mpz_get_ui(r1);
		modular_roots[i].root2 = mpz_get_ui(r2);
	}
	mpz_clear(tmp);
	mpz_clear(r1);
	mpz_clear(r2);
	STOP_TIMER_PRINT_TIME("\nModular roots calculation done");

	/*for(i=0; i<nb_qr_primes; i++)
	 {
	 printf("[%10" PRId64 "-> roots: %10u - %10u]\n", primes[i], modular_roots[i].root1, modular_roots[i].root2);
	 }*/

//--------------------------------------------------------
//         ***** initialize the matrix *****
//--------------------------------------------------------
	START_TIMER();
	init_matrix(&matrix, nb_qr_primes + NB_VECTORS_OFFSET, nb_qr_primes);
	mpz_init2(tmp_matrix_row, nb_qr_primes);
	STOP_TIMER_PRINT_TIME("\nMatrix initialized");

//--------------------------------------------------------
// [Sieving]
//--------------------------------------------------------
	START_TIMER();

	mpz_t x, sieving_index, next_sieving_index;
	unsigned long ui_index, SIEVING_STEP = 50000; /* we sieve for 50000 elements at each loop */
	uint64_t p_pow;
	smooth_number_t *x_squared;

	x_squared = calloc(SIEVING_STEP, sizeof(smooth_number_t));
	smooth_numbers = calloc(nb_qr_primes + NB_VECTORS_OFFSET,
			sizeof(smooth_number_t));

	mpz_init_set(x, sqrtN);
	mpz_init_set(sieving_index, x);
	mpz_init_set(next_sieving_index, x);

	mpz_t p;
	mpz_init(p);
	mpz_t str;
	mpz_init_set(str, sieving_index);
	printf("\nSieving ...\n");

//--------------------------------------------------------
// Init before sieving
//--------------------------------------------------------
	for (i = 0; i < SIEVING_STEP; i++) {
		mpz_init(x_squared[i].value_x);
		mpz_init(x_squared[i].value_x_squared);

		/* the factors_exp array is used to keep track of exponents */
		//x_squared[i].factors_exp = calloc(nb_qr_primes, sizeof(uint64_t));
		/* we use directly the exponents vector modulo 2 to preserve space */mpz_init2(
				x_squared[i].factors_vect, nb_qr_primes);
		mpz_add_ui(x, x, 1);
	}

	int nb_smooth_per_round = 0;
	char s[512];

//--------------------------------------------------------
// WHILE smooth numbers found less than the primes in the smooth base + NB_VECTORS_OFFSET
//--------------------------------------------------------
	while (nb_smooth_numbers_found < nb_qr_primes + NB_VECTORS_OFFSET) {
		nb_smooth_per_round = 0;
		mpz_set(x, next_sieving_index); /* sieve numbers from sieving_index to sieving_index + sieving_step */
		mpz_set(sieving_index, next_sieving_index);

		printf("\r");
		printf(
				"\t\tSieving at: %s30 <--> Smooth numbers found: %" PRId64 "/%" PRId64 "",
				mpz_get_str(NULL, 10, sieving_index), nb_smooth_numbers_found,
				nb_qr_primes);
		fflush(stdout);

		for (i = 0; i < SIEVING_STEP; i++) {
			mpz_set(x_squared[i].value_x, x);

			mpz_pow_ui(x_squared[i].value_x_squared, x, 2); /* calculate value_x_squared <- x²-n */
			mpz_sub(x_squared[i].value_x_squared, x_squared[i].value_x_squared,
					N);

			mpz_clear(x_squared[i].factors_vect);
			mpz_init2(x_squared[i].factors_vect, nb_qr_primes); /* reconstruct a new fresh 0ed vector of size nb_qr_primes bits */

			mpz_add_ui(x, x, 1);
		}
		mpz_set(next_sieving_index, x);

//--------------------------------------------------------
// eliminate factors in the x_squared array, those who are 'destructed' to 1 are smooth
//--------------------------------------------------------

		for (i = 0; i < nb_qr_primes; i++) {
			mpz_set_ui(p, (unsigned long) primes[i]);
			mpz_set(x, sieving_index);

			/* get the first multiple of p that is directly larger that sieving_index
			 * Quadratic SIEVING: all elements from this number and in positions multiples of root1 and root2
			 * are also multiples of p */
			get_sieving_start_index(x, x, p, modular_roots[i].root1);
			mpz_set(str, x);
			mpz_sub(x, x, sieving_index); /* x contains index of first number that is divisible by p */

			for (j = mpz_get_ui(x); j < SIEVING_STEP; j += primes[i]) {
				p_pow = mpz_remove(x_squared[j].value_x_squared,
						x_squared[j].value_x_squared, p); /* eliminate all factors of p */

				if (p_pow & 1) /* mark bit if odd power of p exists in this x_squared[j] */
				{
					mpz_setbit(x_squared[j].factors_vect, i);
				}

				if (mpz_cmp_ui(x_squared[j].value_x_squared, 1) == 0) {
					save_smooth_number(x_squared[j]);
					nb_smooth_per_round++;
				}
				/* sieve next element located p steps from here */
			}

			/* same goes for root2 */
			if (modular_roots[i].root2 == modular_roots[i].root1)
				continue;

			mpz_set(x, sieving_index);

			get_sieving_start_index(x, x, p, modular_roots[i].root2);
			mpz_set(str, x);
			mpz_sub(x, x, sieving_index);

			for (j = mpz_get_ui(x); j < SIEVING_STEP; j += primes[i]) {
				p_pow = mpz_remove(x_squared[j].value_x_squared,
						x_squared[j].value_x_squared, p);

				if (p_pow & 1) {
					mpz_setbit(x_squared[j].factors_vect, i);
				}

				if (mpz_cmp_ui(x_squared[j].value_x_squared, 1) == 0) {
					save_smooth_number(x_squared[j]);
					nb_smooth_per_round++;
				}
			}
		}
		//printf("\tSmooth numbers found %" PRId64 "\n", nb_smooth_numbers_found);
		/*sprintf(s, "[start: %s - end: %s - step: %" PRId64 "] nb_smooth_per_round: %d",
		 mpz_get_str(NULL, 10, sieving_index),
		 mpz_get_str(NULL, 10, next_sieving_index),
		 SIEVING_STEP,
		 nb_smooth_per_round);
		 APPEND_TO_LOG_FILE(s);*/
	}

	STOP_TIMER_PRINT_TIME("\nSieving DONE");

	uint64_t t = 0;

//--------------------------------------------------------
//the matrix ready, start Gauss elimination. The Matrix is filled on the call of save_smooth_number()
//--------------------------------------------------------
	START_TIMER();
	gauss_elimination(&matrix);
	STOP_TIMER_PRINT_TIME("\nGauss elimination done");
	//print_matrix_matrix(&matrix);
	//print_matrix_identity(&matrix);

	uint64_t row_index = nb_qr_primes + NB_VECTORS_OFFSET - 1; /* last row in the matrix */
	int nb_linear_relations = 0;
	mpz_t linear_relation_z, solution_z;
	mpz_init(linear_relation_z);
	mpz_init(solution_z);

	get_matrix_row(linear_relation_z, &matrix, row_index--); /* get the last few rows in the Gauss eliminated matrix*/
	while (mpz_cmp_ui(linear_relation_z, 0) == 0) {
		nb_linear_relations++;
		get_matrix_row(linear_relation_z, &matrix, row_index--);
	}

	printf("\tLinear dependent relations found : %d\n", nb_linear_relations);

//--------------------------------------------------------
// Factor
//--------------------------------------------------------
	//We use the last linear relation to reconstruct our solution
	START_TIMER();
	printf("\nFactorizing..\n");
	mpz_t solution_X, solution_Y;
	mpz_init(solution_X);
	mpz_init(solution_Y);

	/* we start testing from the first linear relation encountered in the matrix */
	for (j = nb_linear_relations; j > 0; j--) {
		printf("Trying %d..\n", nb_linear_relations - j + 1);
		mpz_set_ui(solution_X, 1);
		mpz_set_ui(solution_Y, 1);

		get_identity_row(solution_z, &matrix,
				nb_qr_primes + NB_VECTORS_OFFSET - j + 1);

		for (i = 0; i < nb_qr_primes; i++) {
			if (mpz_tstbit(solution_z, i)) {
				mpz_mul(solution_X, solution_X, smooth_numbers[i].value_x);
				mpz_mod(solution_X, solution_X, N); /* reduce x to modulo N */

				mpz_mul(solution_Y, solution_Y,
						smooth_numbers[i].value_x_squared);
				/*TODO: handling huge stuff here, there is no modulo N like in the solution_X case!
				 * eliminate squares as long as you go*/
			}
		}

		mpz_sqrt(solution_Y, solution_Y);
		mpz_mod(solution_Y, solution_Y, N); /* y = sqrt(MUL(xi²-n)) mod N */

		mpz_sub(solution_X, solution_X, solution_Y);

		mpz_gcd(solution_X, solution_X, N);

		if (mpz_cmp(solution_X, N) != 0 && mpz_cmp_ui(solution_X, 1) != 0) /* factor can be 1 or N, try another relation */
			break;
	}
	mpz_cdiv_q(solution_Y, N, solution_X);

	printf("\n>>>>>>>>>>> FACTORED %s =\n", mpz_get_str(NULL, 10, N));
	printf("\tFactor 1: %s \n\tFactor 2: %s", mpz_get_str(NULL, 10, solution_X),
			mpz_get_str(NULL, 10, solution_Y));

	/*sprintf(s, "\n>>>>>>>>>>> FACTORED %s =\n", mpz_get_str(NULL, 10, N));
	 APPEND_TO_LOG_FILE(s);
	 sprintf(s, "\tFactor 1: %s \n\tFactor 2: %s", mpz_get_str(NULL, 10, solution_X), mpz_get_str(NULL, 10, solution_Y));
	 APPEND_TO_LOG_FILE(s);

	 gettimeofday(&end_global, NULL);
	 timersub(&end_global, &start_global, &elapsed);
	 sprintf(s, "****** TOTAL TIME: %.3f ms\n", elapsed.tv_sec * 1000 + elapsed.tv_usec / (double) 1000);
	 APPEND_TO_LOG_FILE(s);*/

	STOP_TIMER_PRINT_TIME("\nFactorizing done");

	printf("Cleaning memory..\n");

	/********************** clear the x_squared array **********************/
	for (i = 0; i < SIEVING_STEP; i++) {
		mpz_clear(x_squared[i].value_x);
		mpz_clear(x_squared[i].value_x_squared);
		//free(x_squared[i].factors_exp);
		mpz_clear(x_squared[i].factors_vect);
	}
	free(x_squared);
	/********************** clear the x_squared array **********************/

	free(modular_roots);
	/********************** clear the smooth_numbers array **********************/
	for (i = 0; i < nb_qr_primes + NB_VECTORS_OFFSET; i++) {
		mpz_clear(smooth_numbers[i].value_x);
		mpz_clear(smooth_numbers[i].value_x_squared);
		//free(smooth_numbers[i].factors_exp);
	}
	free(smooth_numbers);
	/********************** clear the smooth_numbers array **********************/

	free(primes);
	/********************** clear mpz _t **********************/mpz_clear(B);
	mpz_clear(N);
	sqrtN, rem;
	mpz_clear(x);
	mpz_clear(sieving_index);
	mpz_clear(next_sieving_index);
	mpz_clear(p);
	mpz_clear(str);
	/********************** clear mpz _t **********************/

	free_matrix(&matrix);

	gettimeofday(&end_global, NULL);
	timersub(&end_global, &start_global, &elapsed);
	printf("****** TOTAL TIME: %.3f ms\n",
			elapsed.tv_sec * 1000 + elapsed.tv_usec / (double) 1000);
	show_mem_usage();
	return 0;
}

