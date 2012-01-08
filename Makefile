
all: qs

qs: gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c qs.c
	gcc gmp_patch.h gmp_patch.c eratosthenes.h eratosthenes.c matrix.h matrix.c qs.c -o qs -lgmp -lmpfr -lm -O3

matrix: gmp_patch.h gmp_patch.c matrix.c
	 gcc gmp_patch.h gmp_patch.c matrix.c -o matrix -lgmp -lmpfr -lm
	 
run: qs
	./qs
	
clean:
	rm qs