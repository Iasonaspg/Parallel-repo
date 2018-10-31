#!/bin/bash

for i in *.out
do
    echo "$i results:" >> ../Execution/results/results_$i
    echo "==========================================" >> ../Execution/results/results_$i

    for j in {16..24}
    do
        for z in {1..8}
        do
            printf "%s\t\t" "Input size:" $j >> ../Execution/results/results_$i
            ./$i $j $z 2 >> ../Execution/results/results_$i
            printf "\n" >> ../Execution/results/results_$i
        done
    done
done
