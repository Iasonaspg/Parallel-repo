### Exc1 Structure

Exc1 has 5 folders:

1. **Serial Code**: Contains the source and header files for the sequential Quicksort.
2. **Cilk Code**: Contains the source and header files for the parallel Quicksort using *Intel's Cilk*.
3. **OpenMP Code**: Contains the source and header files for the parallel Quickosrt using *OpenMP*.
4. **Pthreads Code**: Contains the source and header files for the parallel Quicksort using *Pthreads*.
5. **Execution**: Contains the executables for each one parallel and sequential solution, a folder with in-depth execution results and a script that I wrote to produce them.

### Compilation

If you want to recompile each program you should use the appropriate Makefile that I have created. Download the whole repo as zip or clone it and execute the `make` command inside the folder of your interest.

### Results.sh 

*Results.sh* is a shell-script that I created which executes the given programs with an input range of 2^(16-24) array size and 2^(1-8) threads. The produced results are stored in different files that are created by the same script.

If you want to use the same script on your pc follow these steps:

* Create a folder with a name of your choice
* Paste the .out executables that you want to test and the *script.sh*
* Create a new folder in the same directory with the name **results**
* Execute the command `bash script.sh`
* The desired results will be inside the **results** folder

### System Setup

The CPU of the system that the tests took place is an AMD FX-6300 with 6 physical cores and no hyperthreading.

### Custom seed

In C, rand() function produces a pseudo-random sequence of numbers. That sequence is based on the seed, a number that can be defined through the srand() function. If srand() is not used, the default value is srand(1). The third argument of my algorithm is a custom seed. If we use the same seed across all the three programs, the produced array will be the same and the results will be more accurate.
