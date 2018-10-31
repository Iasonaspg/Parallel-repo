
GCC = gcc
CFLAGS = -O3 -pthread

CMAIN=gpthreads.out

all: qsort_pthreads_main.o qsort_pthreads.o
	$(GCC) $(CFLAGS) $^ -o $(CMAIN)

%.o: %.c
	$(GCC) -c $(CFLAGS) $^

clean:
	rm -f *.o *~ $(CMAIN)
