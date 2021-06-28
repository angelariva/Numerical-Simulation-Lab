#!/bin/bash

rm -r results

mkdir -p results

./main 1000 150 square

./main 1000 150 circumference
