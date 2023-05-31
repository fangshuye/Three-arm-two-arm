#!/bin/bash

for r in {4..9}
do
	sbatch Two_arm_indirect.sh $r
done



