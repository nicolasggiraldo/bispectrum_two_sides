CC = mpicc
CFLAGS = -Wall -O3 -I/home/$(USER)/libs/include/
LFLAGS = -L/home/$(USER)/libs/lib -lgsl -lgslcblas -lm -lfftw3 -lconfig
PROGRAM = bispectrum

$(PROGRAM):
	$(CC) -c $@.c $(CFLAGS)
	$(CC) $@.o $(LFLAGS) -o $@.x
	rm $@.o

clean:
	rm -rf *.out
	rm -rf *-
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm *.x
