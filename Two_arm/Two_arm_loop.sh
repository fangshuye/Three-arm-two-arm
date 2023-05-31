#!/bin/bash

for r in {1..9}
do
	sbatch Two_arm.sh $r
done



