/*
 * qs_distributed.c
 *
 *  Created on: Jan 12, 2012
 *      Author: martani
 */
/*-------------- Factorize numbers with the a distributed version of the Quadratic Sieve
 *-------------- Works only with 2 or more processes */

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <inttypes.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <mpi.h>

#include "custom_types.h"
#include "eratosthenes.h"
#include "gmp_patch.h"
#include "matrix.h"
#include "gmp_mpi_lib.h"

struct timeval start_global;
struct timeval end_global;
struct timeval start;
struct timeval end;
struct timeval elapsed;

//only master prints to screen
#define PRINT(id, fmt, ...) \
            do { if (id==0) fprintf(stdout, fmt, __VA_ARGS__); } while (0)

mpz_t N; /* number to factorize */
matrix_t matrix;

uint64_t nb_smooth_numbers_found = 0; /* local for each worker */
uint64_t nb_global_smooth_numbers_found = 0; /* global for everyone*/
uint64_t nb_qr_primes; /* number of primes p where N is a quadratic residue mod p */
uint64_t *primes; /* array holding the primes of the smoothness base */

smooth_number_t *smooth_numbers;
int NB_VECTORS_OFFSET = 5; /* number of additional rows in the matrix, to make sure that a linear relation exists */

smooth_number_t *temp_slaves_smooth_numbers;
int nb_temp_smooth_numbers = 0;

char processor_name[MPI_MAX_PROCESSOR_NAME];

int my_rank, mpi_group_size;
int tag = 0;

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

// save the sooth number n to the smooth_numbers array, and at the same time its exponents vector to the matrix
mpz_t tmp_matrix_row;
void save_smooth_number(smooth_number_t n) {
	if (nb_global_smooth_numbers_found > nb_qr_primes + NB_VECTORS_OFFSET - 1) /* if we have sufficient smooth numbers, skip saving */
		return;

	mpz_clear(tmp_matrix_row); /* tmp_matrix_row must be initialized already */
	mpz_init2(tmp_matrix_row, nb_qr_primes); /* init a vector of *exactly* nb_qr_primes bits */

	smooth_number_t tmp;
	mpz_init(tmp.value_x);
	mpz_init(tmp.value_x_squared);
	mpz_init(tmp.factors_vect);

	mpz_set(tmp.value_x, n.value_x);
	mpz_pow_ui(tmp.value_x_squared, n.value_x, 2);
	mpz_sub(tmp.value_x_squared, tmp.value_x_squared, N);

	mpz_set(tmp.factors_vect, n.factors_vect);

	if (my_rank == 0) /* master appends directly to the matrix */
	{
		smooth_numbers[nb_global_smooth_numbers_found++] = tmp;
		push_row(&matrix, n.factors_vect);
	} else /* append smooth number to the temporary list, later they will be sent to master */
	{
		temp_slaves_smooth_numbers[nb_temp_smooth_numbers++] = tmp;
	}
}

/* called by master to invoke results from slaves */
void gather_smooth_numbers() {
	int i, j;
	for (i = 0; i < mpi_group_size; i++) {
		if (i == my_rank)
			continue;

		temp_slaves_smooth_numbers = gmp_mpi_recv_smooth_numbers(i, tag,
				MPI_COMM_WORLD, &nb_temp_smooth_numbers);

		//PRINT(my_rank, "\n\nreceived %u\n\n", nb_temp_smooth_numbers);

		for (j = 0; j < nb_temp_smooth_numbers; j++) {
			save_smooth_number(temp_slaves_smooth_numbers[j]);
		}
	}
}

/* called by slaves to send their results to master */
void send_smooth_numbers_to_master() {
	gmp_mpi_send_smooth_numbers(temp_slaves_smooth_numbers,
			nb_temp_smooth_numbers, 0, tag, MPI_COMM_WORLD);
	nb_temp_smooth_numbers = 0;
}

/* master sends the number of the global smooth numbers collected,
 * job finishes when enough is collected */
void notify_slaves() {
	unsigned long v = (unsigned long) nb_global_smooth_numbers_found;
	int i;

	for (i = 0; i < mpi_group_size; i++) {
		if (i == my_rank)
			continue;
		MPI_Send(&v, 1, MPI_UNSIGNED_LONG, i, tag, MPI_COMM_WORLD);
	}
}

/* called by slave to get the global number of smooth numbers collected by the group */
unsigned long get_server_notification() {
	unsigned long v;
	MPI_Status status;
	MPI_Recv(&v, 1, MPI_UNSIGNED_LONG, 0, tag, MPI_COMM_WORLD, &status);

	return v;
}

void print_lib_version() {
	PRINT(my_rank, "GMP : %s\n", gmp_version);
	PRINT(my_rank, "MPFR library: %-12s\n\n", mpfr_get_version ());
}

// seems to not work properly if user doesn't have enough privileges
void show_mem_usage() {
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	PRINT(my_rank, "\nUsing: %.2fMB of memory\n", ru.ru_maxrss / 1024.0);
}

void START_TIMER() {
	gettimeofday(&start, NULL);
}

void STOP_TIMER_PRINT_TIME(char *s) {
	gettimeofday(&end, NULL);
	timersub(&end, &start, &elapsed);
	if (my_rank == 0) /* only master prints results to screen */
		PRINT(my_rank, "%s [Time: %.3f ms]\n",
				s, elapsed.tv_sec * 1000 + elapsed.tv_usec / (double) 1000);
}

FILE *fp;
void OPEN_LOG_FILE(char *file_name) {
	char s[512];
	sprintf(s, "%s_%s_distributed", file_name, mpz_get_str(NULL, 10, N));
	if (my_rank != 0) /* only master logs to the file */
	{
		fp = NULL;
		return;
	}
	if ((fp = fopen(s, "w")) == NULL) {
		PRINT(my_rank, "%s", "Cannot open file.\n");
		fp = NULL;
	}
}

void APPEND_TO_LOG_FILE(char *text) {
	if (fp != NULL)
		fprintf(fp, "%s\n", text);
}

/* assuming slaves (workers)) are all homogenous, let them all do the calculations
 regarding primes sieving, calculating the smoothness base and the modular roots */
int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_group_size);
	int len;
	MPI_Get_processor_name(processor_name, &len);

	gettimeofday(&start_global, NULL);
	print_lib_version();

	mpz_init(N);
	mpz_t B;
	mpz_init(B);

	unsigned long int uBase;
	int64_t nb_primes;
	modular_root_t *modular_roots;

	uint64_t i, j;

	if (argc < 2) {
		PRINT(my_rank, "usage: %s Number_to_factorize\n", argv[0]);
		exit(2);
	}

	if (mpz_init_set_str(N, argv[1], 10) == -1) {
		PRINT(my_rank, "Cannot load N %s\n", argv[1]);
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
		PRINT(my_rank, "\n<<<[FACTOR]>>> %s\n", mpz_get_str(NULL, 10, sqrtN));
		return 0;
	}

	if (mpz_probab_prime_p(N, 10) > 0) /* don't bother factoring */
	{
		PRINT(my_rank, "N:%s is prime\n", mpz_get_str(NULL, 10, N));
		exit(0);
	}

	OPEN_LOG_FILE("freq");

//--------------------------------------------------------
//  calculate the smoothness base for the given N
//--------------------------------------------------------
	get_smoothness_base(B, N); /* if N is too small, the program will surely fail, please consider a pen and paper instead */
	uBase = mpz_get_ui(B);
	PRINT(my_rank, "n: %s\tBase: %s\n",
			mpz_get_str(NULL, 10, N), mpz_get_str(NULL, 10, B));

//--------------------------------------------------------
// sieve primes that are less than the smoothness base using Eratosthenes sieve
//--------------------------------------------------------
	START_TIMER();
	nb_primes = sieve_primes_up_to((int64_t) (uBase));

	PRINT(my_rank, "\tPrimes found %" PRId64 " [Smoothness Base %lu]\n",
			nb_primes, uBase);
	STOP_TIMER_PRINT_TIME("\tEratosthenes Sieving done");

//--------------------------------------------------------
// fill the primes array with primes to which n is a quadratic residue
//--------------------------------------------------------
	START_TIMER();
	primes = calloc(nb_primes, sizeof(int64_t));
	nb_qr_primes = fill_primes_with_quadratic_residue(primes, N);

	/*for(i=0; i<nb_qr_primes; i++)
	 PRINT(my_rank, "%" PRId64 "\n", primes[i]);*/

	PRINT(my_rank, "\tN-Quadratic primes found %" PRId64 "\n", nb_qr_primes);
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
	STOP_TIMER_PRINT_TIME("Modular roots calculation done");

//--------------------------------------------------------
//         ***** initialize the matrix *****
//--------------------------------------------------------
	if (my_rank == 0) /* only the master have the matrix */
	{
		START_TIMER();
		init_matrix(&matrix, nb_qr_primes + NB_VECTORS_OFFSET, nb_qr_primes);
		mpz_init2(tmp_matrix_row, nb_qr_primes);
		STOP_TIMER_PRINT_TIME("Matrix initialized");
	}

//--------------------------------------------------------
// [Sieving] - everyones sieves including the master
//--------------------------------------------------------
	START_TIMER();

	mpz_t x, sieving_index, next_sieving_index, relative_start, global_step;
	unsigned long ui_index, SIEVING_STEP = 50000; /* we sieve for 50000 elements at each loop */
	int LOCAL_SIEVING_ROUNDS = 10; /* number of iterations a worker sieves before communicating results to the master */
	unsigned long sieving_round = 0;
	unsigned long nb_big_rounds = 0;

	uint64_t p_pow;
	smooth_number_t *x_squared;

	x_squared = calloc(SIEVING_STEP, sizeof(smooth_number_t));

	if (my_rank == 0)
		smooth_numbers = calloc(nb_qr_primes + NB_VECTORS_OFFSET,
				sizeof(smooth_number_t));
	else
		temp_slaves_smooth_numbers = calloc(500, sizeof(smooth_number_t));
	/* TODO: this is not properly correct, using a linkedlist is better to keep track of temporary
	 * smooth numbers at the slaves nodes however it's pretty rare to find 500 smooth numbers in
	 * 50000 * 10 interval. */

	mpz_init_set(x, sqrtN);
	mpz_init(global_step);
	mpz_init(relative_start);
	mpz_init(sieving_index);
	mpz_init(next_sieving_index);

	mpz_t p;
	mpz_init(p);
	mpz_t str;
	mpz_init_set(str, sieving_index);
	PRINT(my_rank, "\n[%s] Sieving ...\n", processor_name);

//--------------------------------------------------------
// Init before sieving
//--------------------------------------------------------
	for (i = 0; i < SIEVING_STEP; i++) {
		mpz_init(x_squared[i].value_x);
		mpz_init(x_squared[i].value_x_squared);

		mpz_init2(x_squared[i].factors_vect, nb_qr_primes);
		mpz_add_ui(x, x, 1);
	}

	int nb_smooth_per_round = 0;
	char s[512];

//--------------------------------------------------------
// WHILE smooth numbers found less than the primes in the smooth base + NB_VECTORS_OFFSET for master
// Or master asked for more smooth numbers from slaves
//--------------------------------------------------------
	while (1) {
		mpz_set_ui(global_step, nb_big_rounds); /* calculates the coordinate where the workers start sieving from */
		mpz_mul_ui(global_step, global_step, (unsigned long) mpi_group_size);
		mpz_mul_ui(global_step, global_step, SIEVING_STEP);
		mpz_mul_ui(global_step, global_step, LOCAL_SIEVING_ROUNDS);
		mpz_add(global_step, global_step, sqrtN);

		mpz_set_ui(relative_start, SIEVING_STEP);
		mpz_mul_ui(relative_start, relative_start, LOCAL_SIEVING_ROUNDS);
		mpz_mul_ui(relative_start, relative_start, (unsigned long) my_rank);
		mpz_add(relative_start, relative_start, global_step);

		mpz_set(sieving_index, relative_start);
		mpz_set(next_sieving_index, relative_start);

		for (sieving_round = 0; sieving_round < LOCAL_SIEVING_ROUNDS; /* each slave sieves for LOCAL_SIEVING_ROUNDS rounds */
		sieving_round++) {
			nb_smooth_per_round = 0;
			mpz_set(x, next_sieving_index); /* sieve numbers from sieving_index to sieving_index + sieving_step */
			mpz_set(sieving_index, next_sieving_index);

			if (my_rank == 0) {
				printf("\r");
				printf(
						"\t\tSieving at: %s30 <--> Smooth numbers found: %" PRId64 "/%" PRId64 "",
						mpz_get_str(NULL, 10, sieving_index),
						nb_global_smooth_numbers_found, nb_qr_primes);
				fflush(stdout);
			}

			for (i = 0; i < SIEVING_STEP; i++) {
				mpz_set(x_squared[i].value_x, x);

				mpz_pow_ui(x_squared[i].value_x_squared, x, 2); /* calculate value_x_squared <- x²-n */
				mpz_sub(x_squared[i].value_x_squared,
						x_squared[i].value_x_squared, N);

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
		}

		if (my_rank == 0) /* master gathers smooth numbers from slaves */
		{
			gather_smooth_numbers();
			notify_slaves();
		} else /* slaves send their smooth numbers to master */
		{
			send_smooth_numbers_to_master();
			nb_global_smooth_numbers_found = get_server_notification();
		}

		if (nb_global_smooth_numbers_found >= nb_qr_primes + NB_VECTORS_OFFSET)
			break;

		nb_big_rounds++;
	}

	STOP_TIMER_PRINT_TIME("\nSieving DONE");

	if (my_rank == 0) {
		uint64_t t = 0;

//--------------------------------------------------------
//the matrix ready, start Gauss elimination. The Matrix is filled on the call of save_smooth_number()
//--------------------------------------------------------
		START_TIMER();
		gauss_elimination(&matrix);
		STOP_TIMER_PRINT_TIME("\nGauss elimination done");

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

		PRINT(my_rank, "\tLinear dependent relations found : %d\n",
				nb_linear_relations);

//--------------------------------------------------------
// Factor
//--------------------------------------------------------
		//We use the last linear relation to reconstruct our solution
		START_TIMER();
		PRINT(my_rank, "%s", "\nFactorizing..\n");
		mpz_t solution_X, solution_Y;
		mpz_init(solution_X);
		mpz_init(solution_Y);

		/* we start testing from the first linear relation encountered in the matrix */
		for (j = nb_linear_relations; j > 0; j--) {
			PRINT(my_rank, "Trying %d..\n", nb_linear_relations - j + 1);
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

		PRINT(my_rank, "\n>>>>>>>>>>> FACTORED %s =\n",
				mpz_get_str(NULL, 10, N));
		PRINT(
				my_rank,
				"\tFactor 1: %s \n\tFactor 2: %s",
				mpz_get_str(NULL, 10, solution_X), mpz_get_str(NULL, 10, solution_Y));

		sprintf(s, "\n>>>>>>>>>>> FACTORED %s =\n", mpz_get_str(NULL, 10, N));
		APPEND_TO_LOG_FILE(s);
		sprintf(s, "\tFactor 1: %s \n\tFactor 2: %s",
				mpz_get_str(NULL, 10, solution_X),
				mpz_get_str(NULL, 10, solution_Y));
		APPEND_TO_LOG_FILE(s);

		gettimeofday(&end_global, NULL);
		timersub(&end_global, &start_global, &elapsed);
		sprintf(s, "****** TOTAL TIME: %.3f ms\n",
				elapsed.tv_sec * 1000 + elapsed.tv_usec / (double) 1000);
		APPEND_TO_LOG_FILE(s);

		STOP_TIMER_PRINT_TIME("\nFactorizing done");
	}

	PRINT(my_rank, "%s", "\nCleaning memory..\n");

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
	if (my_rank == 0) {
		for (i = 0; i < nb_qr_primes + NB_VECTORS_OFFSET; i++) {
			mpz_clear(smooth_numbers[i].value_x);
			mpz_clear(smooth_numbers[i].value_x_squared);
			mpz_clear(smooth_numbers[i].factors_vect);
			//free(smooth_numbers[i].factors_exp);
		}
		free(smooth_numbers);
	} else {
		for (i = 0; i < 500; i++) {
			mpz_clear(temp_slaves_smooth_numbers[i].value_x);
			mpz_clear(temp_slaves_smooth_numbers[i].value_x_squared);
			mpz_clear(temp_slaves_smooth_numbers[i].factors_vect);
		}
		free(temp_slaves_smooth_numbers);
	}
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
	PRINT(my_rank, "****** TOTAL TIME: %.3f ms\n",
			elapsed.tv_sec * 1000 + elapsed.tv_usec / (double) 1000);
	show_mem_usage();

	MPI_Finalize();

	return 0;
}
