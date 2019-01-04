#!/bin/bash
#PBS -N easyDel
#PBS -j oe
#PBS -q pdlab 
#PBS -l nodes=2:ppn=4,walltime=3:00
cd $PBS_O_WORKDIR

module load gcc/7.2.0
module load openmpi/3.0.0

mpicc MpiFindMedian.c -o mpi -Wall -g
mpiexec -np $PBS_NP mpi 20  > results-8.txt
