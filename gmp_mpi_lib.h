/*
 * gmp_mpi_lib.h
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

#include "custom_types.h"

/* sends an array of type smooth_number_t using MPI */
void gmp_mpi_send_smooth_numbers(smooth_number_t *, int nb_elts, int dest, int tag, MPI_Comm comm);

/* reconstructs received smooth_number_t array, returns how many elements were received in len */
smooth_number_t* gmp_mpi_recv_smooth_numbers(int src, int tag, MPI_Comm comm, int *len);

