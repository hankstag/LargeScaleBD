#! /bin/bash

for i in $( ls LSCM_models )
do 
	echo "deal with " $i
	sbatch --mem=10000 --time=1:00:00 job.s $i
done
