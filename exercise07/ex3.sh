#!/bin/bash

rm -r ex7.4
mkdir -p ex7.4

cp code_MD/input.dummy code_MD/input.gas
sed -i 's/temp/1.2/g' code_MD/input.gas
sed -i 's/npart/108/g' code_MD/input.gas
sed -i 's/rho/0.05/g' code_MD/input.gas
sed -i 's/rcut/5.0/g' code_MD/input.gas
sed -i 's/delta/0.0005/g' code_MD/input.gas
sed -i 's/nstep/1000/g' code_MD/input.gas
sed -i 's/iprint/1000/g' code_MD/input.gas
sed -i 's/measure_time_interval/10/g' code_MD/input.gas
sed -i 's/n_blocks/1/g' code_MD/input.gas

cp code_MD/input.dummy code_MD/input.liquid
sed -i 's/temp/1.1/g' code_MD/input.liquid
sed -i 's/npart/108/g' code_MD/input.liquid
sed -i 's/rho/0.8/g' code_MD/input.liquid
sed -i 's/rcut/2.5/g' code_MD/input.liquid
sed -i 's/delta/0.0005/g' code_MD/input.liquid
sed -i 's/nstep/1000/g' code_MD/input.liquid
sed -i 's/iprint/1000/g' code_MD/input.liquid
sed -i 's/measure_time_interval/10/g' code_MD/input.liquid
sed -i 's/n_blocks/1/g' code_MD/input.liquid

cp code_MD/input.dummy code_MD/input.solid
sed -i 's/temp/0.8/g' code_MD/input.solid
sed -i 's/npart/108/g' code_MD/input.solid
sed -i 's/rho/1.1/g' code_MD/input.solid
sed -i 's/rcut/2.2/g' code_MD/input.solid
sed -i 's/delta/0.0005/g' code_MD/input.solid
sed -i 's/nstep/1000/g' code_MD/input.solid
sed -i 's/iprint/1000/g' code_MD/input.solid
sed -i 's/measure_time_interval/10/g' code_MD/input.solid
sed -i 's/n_blocks/1/g' code_MD/input.solid

for state in solid liquid gas
do
  mkdir -p ex7.4/$state
  mkdir -p ex7.4/$state/results
  mkdir -p ex7.4/$state/frames
  cp code_MD/Primes ex7.4/$state/Primes
  cp code_MD/seed.in ex7.4/$state/seed.in
  cp code_MD/input.$state ex7.4/$state/input.$state
  cp code_MD/config.0 ex7.4/$state/config.0
  cp code_MD/main ex7.4/$state
  cp code_MD/clean.sh ex7.4/$state
  cd ex7.4/$state
  for ((irnd=1; irnd<7; ++irnd))
  do
    if [[ "$irnd" -eq "1" ]]
    then
      ./main --equilibration input.$state
    else
      ./main --restart input.$state
    fi
    mv results results$irnd
    mkdir -p results
  done
  cd ..
  cd ..
done

cp code_MD/input.dummy ex7.4/gas/input.gas
sed -i 's/temp/1.2/g' ex7.4/gas/input.gas
sed -i 's/npart/108/g' ex7.4/gas/input.gas
sed -i 's/rho/0.05/g' ex7.4/gas/input.gas
sed -i 's/rcut/5.0/g' ex7.4/gas/input.gas
sed -i 's/delta/0.0005/g' ex7.4/gas/input.gas
sed -i 's/nstep/40000/g' ex7.4/gas/input.gas
sed -i 's/iprint/1000/g' ex7.4/gas/input.gas
sed -i 's/measure_time_interval/1/g' ex7.4/gas/input.gas
sed -i 's/n_blocks/20/g' ex7.4/gas/input.gas

cp code_MD/input.dummy ex7.4/liquid/input.liquid
sed -i 's/temp/1.1/g' ex7.4/liquid/input.liquid
sed -i 's/npart/108/g' ex7.4/liquid/input.liquid
sed -i 's/rho/0.8/g' ex7.4/liquid/input.liquid
sed -i 's/rcut/2.5/g' ex7.4/liquid/input.liquid
sed -i 's/delta/0.0005/g' ex7.4/liquid/input.liquid
sed -i 's/nstep/40000/g' ex7.4/liquid/input.liquid
sed -i 's/iprint/1000/g' ex7.4/liquid/input.liquid
sed -i 's/measure_time_interval/1/g' ex7.4/liquid/input.liquid
sed -i 's/n_blocks/20/g' ex7.4/liquid/input.liquid

cp code_MD/input.dummy ex7.4/solid/input.solid
sed -i 's/temp/0.8/g' ex7.4/solid/input.solid
sed -i 's/npart/108/g' ex7.4/solid/input.solid
sed -i 's/rho/1.1/g' ex7.4/solid/input.solid
sed -i 's/rcut/2.2/g' ex7.4/solid/input.solid
sed -i 's/delta/0.0005/g' ex7.4/solid/input.solid
sed -i 's/nstep/40000/g' ex7.4/solid/input.solid
sed -i 's/iprint/1000/g' ex7.4/solid/input.solid
sed -i 's/measure_time_interval/1/g' ex7.4/solid/input.solid
sed -i 's/n_blocks/20/g' ex7.4/solid/input.solid

for state in solid liquid gas
do
  cd ex7.4/$state
  ./main --restart input.$state
  cd ..
  cd ..
done
