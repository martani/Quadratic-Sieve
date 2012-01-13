/*
 * gmp_mpi_lib.c
 *
 *  Created on: Jan 12, 2012
 *      Author: martani
 */

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <inttypes.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <mpi.h>
#include <string.h>

#include "gmp_mpi_lib.h"

/* sends an array of type smooth_number_t using MPI
 * <get some popcorn while figuring out how this works, I won't blame you if you give up>
 */
void gmp_mpi_send_smooth_numbers(smooth_number_t *numbers, int nb_elts,
		int dest, int tag, MPI_Comm comm) {
	/*1. send how many elements in the array
	 2. send how many bytes (chars) in the data that was serialized
	 3. send the lengths as a tuples (len(value_x), len(value_x_squared), len(factors_vect))
	 4. send actual data.
	 */
	unsigned long total_buffer_size = 0;
	int sizes[nb_elts * 3];

	char *buffer, *tmp_buffer;
	char *s;
	int t1, t2, t3;
	int i;

	for (i = 0; i < nb_elts; i++) {
		t1 = mpz_sizeinbase(numbers[i].value_x, 10);
		t2 = mpz_sizeinbase(numbers[i].value_x_squared, 10);
		t3 = mpz_sizeinbase(numbers[i].factors_vect, 10);
		sizes[i * 3 + 0] = t1;
		sizes[i * 3 + 1] = t2;
		sizes[i * 3 + 2] = t3;

		total_buffer_size += t1 + t2 + t3;
	}

	tmp_buffer = buffer = malloc(total_buffer_size * sizeof(char) + 1);

	for (i = 0; i < nb_elts; i++) {
		s = mpz_get_str(NULL, 10, numbers[i].value_x);
		strcpy(tmp_buffer, s);
		if (s[sizes[i * 3 + 0] - 1] == '\0') /* gmp adds some noise in the case of some numbers; */
		{
			sizes[i * 3 + 0]--;
			total_buffer_size--;
		}
		tmp_buffer += sizes[i * 3 + 0]; /* advance to the end of the just copied string */
		free(s);

		s = mpz_get_str(NULL, 10, numbers[i].value_x_squared);
		strcpy(tmp_buffer, s);
		if (s[sizes[i * 3 + 1] - 1] == '\0') {
			sizes[i * 3 + 1]--;
			total_buffer_size--;
		}
		tmp_buffer += sizes[i * 3 + 1];
		free(s);

		s = mpz_get_str(NULL, 10, numbers[i].factors_vect);
		strcpy(tmp_buffer, s);
		if (s[sizes[i * 3 + 2] - 1] == '\0') {
			sizes[i * 3 + 2]--;
			total_buffer_size--;
		}
		tmp_buffer += sizes[i * 3 + 2];
		free(s);
	}
	buffer[total_buffer_size] = '\0';

	total_buffer_size++;
	MPI_Send(&nb_elts, 1, MPI_INT, dest, tag, comm);
	MPI_Send(&total_buffer_size, 1, MPI_UNSIGNED_LONG, dest, tag, comm);
	MPI_Send(sizes, nb_elts * 3, MPI_INT, dest, tag, comm);
	MPI_Send(buffer, total_buffer_size, MPI_CHAR, dest, tag, comm);

	free(buffer);
}

/* reconstructs received smooth_number_t elements, returns how many of them was received */
smooth_number_t* gmp_mpi_recv_smooth_numbers(int src, int tag, MPI_Comm comm,
		int *len) {
	MPI_Status status;
	int nb_elts;
	unsigned long total_buffer_size;
	int *sizes;
	char *buffer, *old_buffer;
	char *s;

	int i;

	/* get length information */
	MPI_Recv(&nb_elts, 1, MPI_INT, src, tag, comm, &status);
	MPI_Recv(&total_buffer_size, 1, MPI_UNSIGNED_LONG, src, tag, comm, &status);

	sizes = malloc(nb_elts * 3 * sizeof(int));
	old_buffer = buffer = malloc(total_buffer_size * sizeof(char));

	smooth_number_t *numbers = malloc(nb_elts * sizeof(smooth_number_t));

	MPI_Recv(sizes, nb_elts * 3, MPI_INT, src, tag, comm, &status);
	MPI_Recv(buffer, total_buffer_size, MPI_CHAR, src, tag, comm, &status);

	for (i = 0; i < nb_elts; i++) {
		mpz_init(numbers[i].value_x);
		mpz_init(numbers[i].value_x_squared);
		mpz_init(numbers[i].factors_vect);

		s = malloc(sizes[i * 3 + 0] + 1);
		strncpy(s, buffer, sizes[i * 3 + 0]);
		s[sizes[i * 3 + 0]] = '\0';
		mpz_set_str(numbers[i].value_x, s, 10);
		free(s);
		buffer += sizes[i * 3 + 0];

		s = malloc(sizes[i * 3 + 1] + 1);
		strncpy(s, buffer, sizes[i * 3 + 1]);
		s[sizes[i * 3 + 1]] = '\0';
		mpz_set_str(numbers[i].value_x_squared, s, 10);
		free(s);
		buffer += sizes[i * 3 + 1];

		s = malloc(sizes[i * 3 + 2] + 1);
		strncpy(s, buffer, sizes[i * 3 + 2]);
		s[sizes[i * 3 + 2]] = '\0';
		mpz_set_str(numbers[i].factors_vect, s, 10);
		free(s);
		buffer += sizes[i * 3 + 2];
	}
	free(sizes);
	free(old_buffer);

	*len = nb_elts;
	return numbers;
}

