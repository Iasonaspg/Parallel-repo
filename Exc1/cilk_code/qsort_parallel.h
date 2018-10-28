/*

Parallel implemantation of quicksort using Intel Cilk Plus

Author: Michail Iason Pavlidis
Email: michailpg@ece.auth.gr

*/

#ifndef _QSORT_PARALLEL_H_
#define _QSORT_PARALLEL_H_

/* Sort vector v (of length n) using QuickSort */
void qsort_par(int *v, int n);

#endif /* _QSORT-PARALLEL_H_ */
