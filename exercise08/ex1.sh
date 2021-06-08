#!/bin/bash

cd ex8.1

MU_MIN=0.7
MU_MAX=0.9
MU_STEP=0.01

SIGMA_MIN=0.5
SIGMA_MAX=0.8
SIGMA_STEP=0.01

echo "Script to run Variational Monte Carlo"

echo "Generating parameters..."

if [ -f variational.dat ]; then
    rm variational.dat
fi

touch variational.dat

for mu in `seq $MU_MIN $MU_STEP $MU_MAX`
do
  for sigma in `seq $SIGMA_MIN $SIGMA_STEP $SIGMA_MAX`
  do
    echo "$mu $sigma"
    echo ${mu/,/.}
    echo ${sigma/,/.}
    ./main  ${mu/,/.} ${sigma/,/.}
    tail -1 energy.dat | awk -v mu="${mu/,/.}" -v sigma="${sigma/,/.}" '{print mu, sigma, $2, $3}' >> variational.dat
  done
done
# sort --key 3 --numeric-sort variational.dat
RESULTS=($(awk '{print $3, $1, $2}' variational.dat | sort --numeric-sort | head -1))

echo "Optimal parameters are: mu=${RESULTS[1]}, sigma=${RESULTS[2]}"

echo "Results saved in 'variational.dat'..."

echo "All done!"

echo "Now we run a simulation with optimal parameters and 10 000 000 steps"

./main ${RESULTS[1]} ${RESULTS[2]} 10000000

cd ..
