#!/bin/bash

rm -r ex6.1
mkdir -p ex6.1


for state in metropolis gibbs
do
  for h in 0 0.02
  do
    mkdir -p ex6.1/$state
    mkdir -p ex6.1/$state/$h
    mkdir -p ex6.1/$state/${h}/results

    cp code/ex.py ex6.1/ex.py
    cp code/Primes ex6.1/$state/${h}/Primes
    cp code/seed.in ex6.1/$state/${h}/seed.in
    cp code/main ex6.1/$state/${h}/main

  done
done

cd ex6.1
python3 ex.py
