#!/bin/bash

for j in *.out
do
    echo "$j results:" >> results_$j
    echo "==========================================" >> results_$j
    for i in {12..24}
    do
        for x in {0..8}
        do
            if [ "$i" == "24" ]
            then
                end=2
            else 
                end=0
            fi
            for z in $(seq 0 $end)
            do
                printf "%s\t\t" "Input size:" $i >> results_$j
                ./$j $i $x 2 >> results_$j
                printf "\n" >> results_$j
            done
        done
     done
done

echo "Sequential results:" >> results.txt
echo "==========================================" >> results.txt

for j in {12..24}
do
    printf "%s\t\t" "Input size:" $j >> results.txt
    ./qsort $j 2 >> results.txt
    printf "\n" >> results.txt
done
