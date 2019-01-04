#!/bin/bash

for j in {1..1000}
do
    for i in {1..1000}
    do
            seq 60 0.001 65 | shuf -n 1 >> rand4.out
    done
done
