/*

Parallel implemantation of quicksort using POSIX Threads

Author: Michail Iason Pavlidis
Email: michailpg@ece.auth.gr

*/

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "qsort_pthreads.h"


int lvl = 0, f = 1, threads;
pthread_mutex_t mutexsum, mutexsum1;


/* swap -- swap elements k and l of vector v */
void swap(int *v, int k, int l) {
  int temp = v[k];
  v[k] = v[l];
  v[l] = temp;
}


/* partition -- in-place update of elements */
int partition(int *v, int n) {
  int pivot = v[n-1];
  int i = 0;
  
  for (int j = 0; j < n - 1; j++) 
    if (v[j] < pivot) 
      swap(v,i++,j);

  swap(v, i, n - 1);
  return (i);
}

/* Quicksort that is called by threads */
void *qsort_par(void *data) {

  thread_data *td = (thread_data *) data;
  if ( f ){
    pthread_mutex_init(&mutexsum, NULL);
    pthread_mutex_init(&mutexsum1, NULL);
    threads = td->thread_num;
  }

  if ((td->length) > 1) {

    
    int p = partition(td->index, td->length);
    thread_data td1;
    td1.index = &(td->index[p + 1]);
    td1.length = td->length - p -1;
    td->length = p;
           
    
    if ((lvl < threads) && (td->length > 50000)){
      pthread_t thread2;
      pthread_mutex_lock(&mutexsum);
      ++lvl;
      pthread_mutex_unlock(&mutexsum);
      f = 0;
      
      pthread_create(&thread2, NULL, (void *) qsort_par, (void *) td);
      qsort_par(&td1);
      pthread_join(thread2, NULL);

      pthread_mutex_lock(&mutexsum1);
      --lvl;
      pthread_mutex_unlock(&mutexsum1);
      
    }
    else{
      qsort_par(&td1);
      qsort_par(td);
    }
  }
}