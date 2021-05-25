#!/bin/bash

rm -r ex7.2
mkdir -p ex7.2

cp code_MC/input.dummy code_MC/input.gas
sed -i 's/temp/1.2/g' code_MC/input.gas
sed -i 's/npart/108/g' code_MC/input.gas
sed -i 's/rho/0.05/g' code_MC/input.gas
sed -i 's/rcut/5.0/g' code_MC/input.gas
sed -i 's/delta/7.5/g' code_MC/input.gas
sed -i 's/nstep/3000/g' code_MC/input.gas
sed -i 's/n_blocks/1/g' code_MC/input.gas

cp code_MC/input.dummy code_MC/input.liquid
sed -i 's/temp/1.1/g' code_MC/input.liquid
sed -i 's/npart/108/g' code_MC/input.liquid
sed -i 's/rho/0.8/g' code_MC/input.liquid
sed -i 's/rcut/2.5/g' code_MC/input.liquid
sed -i 's/delta/0.2/g' code_MC/input.liquid
sed -i 's/nstep/3000/g' code_MC/input.liquid
sed -i 's/n_blocks/1/g' code_MC/input.liquid

cp code_MC/input.dummy code_MC/input.solid
sed -i 's/temp/0.8/g' code_MC/input.solid
sed -i 's/npart/108/g' code_MC/input.solid
sed -i 's/rho/1.1/g' code_MC/input.solid
sed -i 's/rcut/2.2/g' code_MC/input.solid
sed -i 's/delta/0.12/g' code_MC/input.solid
sed -i 's/nstep/3000/g' code_MC/input.solid
sed -i 's/n_blocks/1/g' code_MC/input.solid

for state in solid liquid gas
do
  mkdir -p ex7.2/$state
  mkdir -p ex7.2/$state/results
  mkdir -p ex7.2/$state/frames
  cp code_MC/Primes ex7.2/$state/Primes
  cp code_MC/seed.in ex7.2/$state/seed.in
  cp code_MC/input.$state ex7.2/$state/input.$state
  cp code_MC/config.0 ex7.2/$state/config.0
  cp code_MC/main ex7.2/$state
  cd ex7.2/$state

  ./main --equilibration input.$state --instant
  cp -r results results_eq
  cd ..
  cd ..
done

cp code_MC/input.dummy code_MC/input.gas
sed -i 's/temp/1.2/g' code_MC/input.gas
sed -i 's/npart/108/g' code_MC/input.gas
sed -i 's/rho/0.05/g' code_MC/input.gas
sed -i 's/rcut/5.0/g' code_MC/input.gas
sed -i 's/delta/7.5/g' code_MC/input.gas
sed -i 's/nstep/2000/g' code_MC/input.gas
sed -i 's/n_blocks/20/g' code_MC/input.gas

cp code_MC/input.dummy code_MC/input.liquid
sed -i 's/temp/1.1/g' code_MC/input.liquid
sed -i 's/npart/108/g' code_MC/input.liquid
sed -i 's/rho/0.8/g' code_MC/input.liquid
sed -i 's/rcut/2.5/g' code_MC/input.liquid
sed -i 's/delta/0.2/g' code_MC/input.liquid
sed -i 's/nstep/2000/g' code_MC/input.liquid
sed -i 's/n_blocks/20/g' code_MC/input.liquid

cp code_MC/input.dummy code_MC/input.solid
sed -i 's/temp/0.8/g' code_MC/input.solid
sed -i 's/npart/108/g' code_MC/input.solid
sed -i 's/rho/1.1/g' code_MC/input.solid
sed -i 's/rcut/2.2/g' code_MC/input.solid
sed -i 's/delta/0.12/g' code_MC/input.solid
sed -i 's/nstep/2000/g' code_MC/input.solid
sed -i 's/n_blocks/20/g' code_MC/input.solid

for state in solid liquid gas
do
  cd ex7.2/$state
  ./main --restart input.$state --instant
  cd ..
  cd ..
done
