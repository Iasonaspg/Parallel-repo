/*

Parallel implemantation of quicksort using POSIX Threads

Author: Michail Iason Pavlidis
Email: michailpg@ece.auth.gr

*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include "qsort_pthreads.h"
#include <pthread.h>

int  test( int *a, int n);
void init( int *a, int n);
void print(int *a, int n);


int main(int argc, char **argv){

  //Parse input
  if ((argc != 4) && (argc != 3)){
      printf("Invalid number of arguments\n");
      exit(1);
  }
  else if ((atoi(argv[1]) == 0) || (atoi(argv[1]) < 12) || (atoi(argv[1]) > 24) || (atoi(argv[2]) > 8)){
      printf("Invalid arguments\n");
      exit(1);
  }


  /* variables to hold execution time */
  struct timeval startwtime, endwtime;
  double par_time;

  /* Memory allocation */
  int n  = 1<<atoi(argv[1]);
  int *a = (int *) malloc(n * sizeof(int));

  if (argc == 4){
      int seed = atoi(argv[3]);
      srand(seed);
  }


  /* initialize vector */
  init(a, n);

  // print(a,n);

  thread_data *tdata;

  tdata = (thread_data *) malloc(sizeof(thread_data));

  tdata->index = a;
  tdata->length = n;
  tdata->thread_num = 1<<atoi(argv[2]);


  /* sort elements */
  gettimeofday(&startwtime, NULL);
  qsort_par(tdata);

  gettimeofday(&endwtime, NULL);

  /* get time in seconds */
  par_time = (double) ((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                    + endwtime.tv_sec - startwtime.tv_sec);

  /* validate result */
  int pass = test(a, n);
  printf(" \tTEST %s\n", (pass) ? "PASSed" : "FAILed");
  //assert( pass != 0 );


  /* print execution time */
  if (atoi(argv[2])){
        printf("\tParallel wall clock time using %d threads: %f sec\n", tdata->thread_num, par_time);
    }
    else{
        printf("\tParallel wall clock time using %d thread: %f sec\n", tdata->thread_num, par_time);
    }

  pthread_exit(NULL);
  free(tdata);
  free(a);

  return 0;

}


/** -------------- SUB-PROCEDURES  ----------------- **/

/** procedure test() : verify sort results **/
int test(int *a, int n) {

  // TODO: implement
  //int pass = 0;

  for (int i = 0; i < n-1; i++) {
    if (*(a + i +1) < *(a + i)) {
      return 0;
    }
  }

  return 1;

}

/** procedure init() : initialize array "a" with data **/
void init(int *a, int n) {
  int i;
  for (i = 0; i < n; i++) {
    a[i] = rand() % n; // (N - i);
  }
}

/** procedure  print() : print array elements **/
void print(int *a, int n) {
  int i;
  for (i = 0; i < n; i++) {
    printf("%d ", a[i]);
  }
  printf("\n");
}
