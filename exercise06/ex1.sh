#!/bin/bash

mkdir -p ex6.1

for state in metropolis gibbs
do
  mkdir -p ex6.1/$state
  for h in 0 0.02
  do
    mkdir -p ex6.1/$state/$h
    mkdir -p ex6.1/$state/${h}/results
    cp code/input.dummy ex6.1/$state/${h}/input.dat
    sed -i "s/h_dummy/0/g" ex6.1/$state/${h}/input.dat
    sed -i "s/nblk_dummy/1/g" ex6.1/$state/${h}/input.dat
    sed -i "s/nstep_dummy/5/g" ex6.1/$state/${h}/input.dat
    if [[ "$state" -eq "metropolis" ]];
    then
      sed -i "s/metro_dummy/1/g" ex6.1/$state/${h}/input.dat
    elif [[ "$state" -eq "gibbs" ]]; then
      sed -i "s/metro_dummy/0/g" ex6.1/$state/${h}/input.dat
    fi

    for t in 2.0 0.5 0.95
    do
      sed -i "s/temp_dummy/${t}/g" ex6.1/$state/${h}/input.dat

      cp code/Primes ex6.1/$state/${h}/Primes
      cp code/seed.in ex6.1/$state/${h}/seed.in
      cp code/main ex6.1/$state/${h}/main

      cd ex6.1/$state/$h
      ./main --equilibration
      cp -r  results results_eq_$t
    done
    cd ..
    cd ..
    cd ..
  done
done


for state in metropolis gibbs
do
  mkdir -p ex6.1/$state
  for h in 0 0.02
  do
    mkdir -p ex6.1/$state/$h
    mkdir -p ex6.1/$state/${h}/results
    sed -i "s/temp_dummy/${t}/g" ex6.1/$state/${h}/input.dat
    mkdir -p ex6.1/$state/${h}
    cp code/input.dummy ex6.1/$state/${h}/input.dat
    sed -i "s/h_dummy/0/g" ex6.1/$state/${h}/input.dat
    sed -i "s/nblk_dummy/2/g" ex6.1/$state/${h}/input.dat
    sed -i "s/nstep_dummy/10/g" ex6.1/$state/${h}/input.dat
    if [[ "$state" -eq "metropolis" ]];
    then
      sed -i "s/metro_dummy/1/g" ex6.1/$state/${h}/input.dat
    elif [[ "$state" -eq "gibbs" ]]; then
      sed -i "s/metro_dummy/0/g" ex6.1/$state/${h}/input.dat
    fi

    for t in 2.0 0.5 0.95
    do
      sed -i "s/temp_dummy/${t}/g" ex6.1/$state/${h}/input.dat

      cp code/Primes ex6.1/$state/${h}/Primes
      cp code/seed.in ex6.1/$state/${h}/seed.in
      cp code/main ex6.1/$state/${h}/main

      cd ex6.1/$state/${h}
      ./main --restart
      cp -r  results results_$t
    done
    cd ..
    cd ..
    cd ..
  done
done
