#!/bin/bash

seed=77

latticegen="/usr/bin/latticegen"

for d in {12..15} #dimension
do
	for i in {0..9} #instance
	do
		((seed=seed+1))
		${latticegen} -randseed ${seed} q $d 1 20 b > "./cartwheel/q_${d}_1_20_b_${i}"
	done
done
