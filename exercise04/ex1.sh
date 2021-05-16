#!/bin/bash


mkdir -p ex4.1

cp code/input.dummy code/input.gas
sed -i 's/temp/1.2/g' code/input.gas
sed -i 's/npart/108/g' code/input.gas
sed -i 's/rho/0.05/g' code/input.gas
sed -i 's/rcut/5.0/g' code/input.gas
sed -i 's/delta/0.0005/g' code/input.gas
sed -i 's/nstep/1000/g' code/input.gas
sed -i 's/iprint/1000/g' code/input.gas
sed -i 's/measure_time_interval/10/g' code/input.gas
sed -i 's/n_blocks/1/g' code/input.gas

cp code/input.dummy code/input.liquid
sed -i 's/temp/1.1/g' code/input.liquid
sed -i 's/npart/108/g' code/input.liquid
sed -i 's/rho/0.8/g' code/input.liquid
sed -i 's/rcut/2.5/g' code/input.liquid
sed -i 's/delta/0.0005/g' code/input.liquid
sed -i 's/nstep/1000/g' code/input.liquid
sed -i 's/iprint/1000/g' code/input.liquid
sed -i 's/measure_time_interval/10/g' code/input.liquid
sed -i 's/n_blocks/1/g' code/input.liquid

cp code/input.dummy code/input.solid
sed -i 's/temp/0.8/g' code/input.solid
sed -i 's/npart/108/g' code/input.solid
sed -i 's/rho/1.1/g' code/input.solid
sed -i 's/rcut/2.2/g' code/input.solid
sed -i 's/delta/0.0005/g' code/input.solid
sed -i 's/nstep/1000/g' code/input.solid
sed -i 's/iprint/1000/g' code/input.solid
sed -i 's/measure_time_interval/10/g' code/input.solid
sed -i 's/n_blocks/1/g' code/input.solid

for state in solid liquid gas
do
  mkdir -p ex4.1/$state
  mkdir -p ex4.1/$state/results
  mkdir -p ex4.1/$state/frames
  cp code/Primes ex4.1/$state/Primes
  cp code/seed.in ex4.1/$state/seed.in
  cp code/input.$state ex4.1/$state/input.$state
  cp code/config.0 ex4.1/$state/config.0
  cp code/main ex4.1/$state
  cp code/clean.sh ex4.1/$state
  cd ex4.1/$state
  for ((irnd=1; irnd<7; ++irnd))
  do
    if [[ "$irnd" -eq "1" ]]
    then
      ./main --equilibration input.$state
    else
      ./main --restart input.$state
    fi
    cp -r results results$irnd
  done
  cd ..
  cd ..
done
