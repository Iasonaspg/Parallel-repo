/*

Parallel implemantation of quicksort using POSIX Threads

Author: Michail Iason Pavlidis
Email: michailpg@ece.auth.gr

*/

#ifndef _QSORT_PTHREADS_H_
#define _QSORT_PTHREADS_H_

/* Sort vector v (of length n) using QuickSort */
void *qsort_par(void *data);

// Structure with data
typedef struct thread_d{
    int*  index;
    int  length;
    int thread_num;
}thread_data;

#endif /* _QSORT-PARALLEL_H_ */
