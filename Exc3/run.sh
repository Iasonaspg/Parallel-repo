#!/bin/bash

for i in $(seq 20 25)
do
    for j in $(seq 4 6)
    do
        ./knn $i $j > results_$i_$j
    done
done
