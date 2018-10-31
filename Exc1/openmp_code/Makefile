
GCC = gcc
CFLAGS = -O3 -fopenmp

CMAIN=gopenmp.out

all: qsort_openmp_main.o qsort_openmp.o
	$(GCC) $(CFLAGS) $^ -o $(CMAIN)

%.o: %.c
	$(GCC) -c $(CFLAGS) $^

clean:
	rm -f *.o *~ $(CMAIN)
