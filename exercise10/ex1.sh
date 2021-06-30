#!/bin/bash

cd ex10.1

rm -r results

mkdir -p results

./main 500 1500 square

./main 500 100 circumference

cd ..
