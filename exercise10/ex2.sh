#!/bin/bash

cd ex10.2

rm -r results

mkdir -p results

mpirun --use-hwthread-cpus  ./main 100 200 square

mpirun --use-hwthread-cpus  ./main 100 200 circumference

cd ..
