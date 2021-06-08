 cd ex8.3
 
 #PIGS
for tau in 05 1 2 3
do
	cp input/input_t$tau.pigs input.dat
	./qmc1d 	#variational trial wave function

	for res in probability potential kinetic
	do
		cp $res.dat results/$res.var_t$tau.dat
	done;echo;echo

done


for tau in 3 5 6
do
	cp input/input_t$tau.pigs input.dat
	./qmc1d_const 	#constant trial wave function

	for res in probability potential kinetic
	do
		cp $res.dat results/$res.const_t$tau.dat
	done;echo;echo

done


#PIMC
for temp in 0_25 1_25 5 50
do
	cp input/input_$temp.pimc input.dat
	./qmc1d

	for res in probability potential kinetic
	do
		cp $res.dat results/$res.T$temp.dat
	done;echo;echo

done

for res in probability potential kinetic
do
	rm -rf $res.dat
done

cd ..
