#!/bin/bash

for r in {2..3}
do
	sbatch Three_arm.sh $r
done



