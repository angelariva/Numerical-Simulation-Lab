#!/bin/bash

rm -r equilibration
mkdir -p equilibration

for state in metropolis gibbs
do
  for h in 0 0.02
  do
    for t in 2.0 0.5 0.95
    do
      mkdir -p equilibration/$state
      mkdir -p equilibration/$state/$h
      mkdir -p equilibration/$state/${h}/${t}
      mkdir -p equilibration/$state/${h}/${t}/results

      cp code/Primes equilibration/$state/${h}/${t}/Primes
      cp code/seed.in equilibration/$state/${h}/${t}/seed.in
      cp code/main equilibration/$state/${h}/${t}/main
      cp code/input.dummy equilibration/$state/${h}/${t}/input.dat

      sed -i "s/h_dummy/${h}/g" equilibration/$state/${h}/${t}/input.dat
      sed -i "s/nblk_dummy/1/g" equilibration/$state/${h}/${t}/input.dat
      sed -i "s/nstep_dummy/5000/g" equilibration/$state/${h}/${t}/input.dat

      if [ "$state" == "metropolis" ];
      then
        sed -i "s/metro_dummy/1/g" equilibration/$state/${h}/${t}/input.dat
      elif [ "$state" == "gibbs" ];
      then
        sed -i "s/metro_dummy/0/g" equilibration/$state/${h}/${t}/input.dat
      fi

      cd equilibration/$state/${h}/${t}
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

: << 'END'
for state in metropolis gibbs
do
  for h in 0 0.02
  do
    for t in 2.0 0.5 0.95
    do
      mkdir -p equilibration/$state
      mkdir -p equilibration/$state/$h
      mkdir -p equilibration/$state/${h}/${t}
      mkdir -p equilibration/$state/${h}/${t}/results

      cp code/input.dummy equilibration/$state/${h}/${t}/input.dat
      cp code/Primes equilibration/$state/${h}/${t}/Primes
      cp code/seed.in equilibration/$state/${h}/${t}/seed.in
      cp code/main equilibration/$state/${h}/${t}/main

      sed -i "s/temp_dummy/${t}/g" equilibration/$state/${h}/${t}/input.dat
      sed -i "s/h_dummy/${h}/g" equilibration/$state/${h}/${t}/input.dat
      sed -i "s/nblk_dummy/20/g" equilibration/$state/${h}/${t}/input.dat
      sed -i "s/nstep_dummy/50000/g" equilibration/$state/${h}/${t}/input.dat
      if [ "$state" == "metropolis" ];
      then
        sed -i "s/metro_dummy/1/g" equilibration/$state/${h}/${t}/input.dat
      elif [ "$state" == "gibbs" ];
      then
        sed -i "s/metro_dummy/0/g" equilibration/$state/${h}/${t}/input.dat
      fi

      cd equilibration/$state/${h}/${t}
      sed -i "s/temp_dummy/${t}/g" input.dat
     ./main --restart
     cd ..
     cd ..
     cd ..
     cd ..
    done
  done
done
END
