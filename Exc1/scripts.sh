#!/bin/bash

for j in *.out
do
    echo "$j results:" >> results.txt
    echo "==========================================" >> results.txt
    for i in {12..24}
    do
        for x in {0..8}
        do
            for z in {0..2}
            do
                printf "%s\t\t" "Input size:" $i >> results.txt
                ./$j $i $x 2 >> results.txt
                printf "\n" >> results.txt
            done
        done
        #echo "---------------------------------------------" >> results.txt
    done
done

for j in {12..24}
do
    printf "%s\t\t" "Input size:" $i >> results.txt
    ./qsort $j 2
