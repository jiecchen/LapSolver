#!/bin/bash

N=1000
T=0

for I in `seq 1 $N`
do
    TIME=$({ time ./test; } |& grep real | cut -dm -f2 | cut -ds -f1)
    T=$(echo $T+$TIME | bc)
done

echo "$T/$N.0" | bc -l
