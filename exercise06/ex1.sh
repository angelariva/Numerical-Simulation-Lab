#!/bin/bash

rm -r ex6.1
mkdir -p ex6.1

for state in metropolis gibbs
do
  for h in 0 0.02
  do
    for t in 2.0 0.5 0.95
    do
      mkdir -p ex6.1/$state
      mkdir -p ex6.1/$state/$h
      mkdir -p ex6.1/$state/${h}/${t}
      mkdir -p ex6.1/$state/${h}/${t}/results

      cp code/Primes ex6.1/$state/${h}/${t}/Primes
      cp code/seed.in ex6.1/$state/${h}/${t}/seed.in
      cp code/main ex6.1/$state/${h}/${t}/main
      cp code/input.dummy ex6.1/$state/${h}/${t}/input.dat

      sed -i "s/h_dummy/0/g" ex6.1/$state/${h}/${t}/input.dat
      sed -i "s/nblk_dummy/1/g" ex6.1/$state/${h}/${t}/input.dat
      sed -i "s/nstep_dummy/5000/g" ex6.1/$state/${h}/${t}/input.dat

      if [ "$state" == "metropolis" ];
      then
        sed -i "s/metro_dummy/1/g" ex6.1/$state/${h}/${t}/input.dat
      elif [ "$state" == "gibbs" ];
      then
        sed -i "s/metro_dummy/0/g" ex6.1/$state/${h}/${t}/input.dat
      fi

      cd ex6.1/$state/${h}/${t}
      sed -i "s/temp_dummy/${t}/g" input.dat
      ./main --equilibration
      mv  results results_eq
      cd ..
      cd ..
      cd ..
      cd ..
      done
  done
done


for state in metropolis gibbs
do
  for h in 0 0.02
  do
    for t in 2.0 0.5 0.95
    do
      mkdir -p ex6.1/$state
      mkdir -p ex6.1/$state/$h
      mkdir -p ex6.1/$state/${h}/${t}
      mkdir -p ex6.1/$state/${h}/${t}/results

      cp code/input.dummy ex6.1/$state/${h}/${t}/input.dat
      cp code/Primes ex6.1/$state/${h}/${t}/Primes
      cp code/seed.in ex6.1/$state/${h}/${t}/seed.in
      cp code/main ex6.1/$state/${h}/${t}/main

      sed -i "s/temp_dummy/${t}/g" ex6.1/$state/${h}/${t}/input.dat
      sed -i "s/h_dummy/0/g" ex6.1/$state/${h}/${t}/input.dat
      sed -i "s/nblk_dummy/20/g" ex6.1/$state/${h}/${t}/input.dat
      sed -i "s/nstep_dummy/50000/g" ex6.1/$state/${h}/${t}/input.dat
      if [ "$state" == "metropolis" ];
      then
        sed -i "s/metro_dummy/1/g" ex6.1/$state/${h}/${t}/input.dat
      elif [ "$state" == "gibbs" ];
      then
        sed -i "s/metro_dummy/0/g" ex6.1/$state/${h}/${t}/input.dat
      fi

      cd ex6.1/$state/${h}/${t}
      sed -i "s/temp_dummy/${t}/g" input.dat
     ./main --restart
     cd ..
     cd ..
     cd ..
     cd ..
    done
  done
done
