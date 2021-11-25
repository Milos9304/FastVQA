#!/bin/bash

for f in rank_2*; do
    if [ -d "$f" ]; then
    # $f is a directory
       rank=${f#rank_*}
       echo $rank
       ../bin/fastVQA --vqe -e -r $rank -x 0.175 &> "log/rank_$rank.log"
    fi
done
