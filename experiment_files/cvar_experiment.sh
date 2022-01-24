#!/bin/bash

for f in rank_15; do
    for x in 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5 0.525 0.55 0.575 0.6 0.625 0.65 0.675 0.7 0.725 0.75 0.775 0.8 0.825 0.85 0.875 0.9 0.925 0.95 0.975 1; do
        if [ -d "$f" ]; then
            # $f is a directory
            rank=${f#rank_*}
            echo rank=$rank cvar=$x
            ../bin/fastVQA --vqe -e -r $rank --prefix cvar$x -x $x &> "log/cvar_${rank}_$x.log"
        fi
    done
done

#above code is previous. below is fix for dim 25
#rank=25
#../bin/fastVQA --vqe -e -r $rank -x 0.175 &> "log/rank_$rank.log"
