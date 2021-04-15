#!/bin/bash

seed=77

latticegen="/usr/bin/latticegen"

for d in {5..10} #dimension
do
	for i in {0..4} #instance
	do
		((seed=seed+1))
		${latticegen} -randseed ${seed} q $d 1 20 b > "./blackeye/q_${d}_1_20_b_${i}"
	done
done
#$latticegen q 10 1 20 b > "q_10_1_20_b"
#$latticegen q 130 1 20 b > "q_130_1_20_b"

