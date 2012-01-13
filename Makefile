
all: qs

qs: gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c qs.c
	gcc gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c qs.c -o qs -lgmp -lmpfr -lm -O3
	 
center: gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c qs_center_n.c
	gcc gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c qs_center_n.c -o qs_center_n -lgmp -lmpfr -lm -O3

basic: gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c qs_exponent_vector.c
	gcc gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c qs_exponent_vector.c -o qs_exponent_vector -lgmp -lmpfr -lm -O3

distributed :gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c custom_types.h gmp_mpi_lib.h gmp_mpi_lib.c qs_distributed.c
	mpicc gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c custom_types.h gmp_mpi_lib.h gmp_mpi_lib.c qs_distributed.c -o qs_distributed -lgmp -lmpfr -lm -O3 

run: qs
	./qs
	
clean:
	rm qs qs_exponent_vector qs_center_n qs_distributed *~