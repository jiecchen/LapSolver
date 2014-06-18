#!/bin/bash

N=100
T=0

for I in `seq 1 $N`
do
    TIME=$({ time ./$1; } 2>&1 | grep real | cut -dm -f2 | cut -ds -f1)
    T=$(echo $T+$TIME | bc)
done

echo "$T/$N.0" | bc -l
