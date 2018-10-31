
GCC = gcc
CFLAGS = -O3 -fcilkplus

CMAIN=gcilk.out

all: qsort_main.o qsort_parallel.o
	$(GCC) $(CFLAGS) $^ -o $(CMAIN)

%.o: %.c
	$(GCC) -c $(CFLAGS) $^

clean:
	rm -f *.o *~ $(CMAIN)
