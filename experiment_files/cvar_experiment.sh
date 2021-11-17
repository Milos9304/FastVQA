#!/bin/bash

for f in r_0.200; do
    if [ -d "$f" ]; then
    # $f is a directory
       ratio=${f#r_*}
       echo $ratio
       ../bin/fastVQA --vqe -e -r 15 -l -x $ratio -f 1 -m 1050 &> "log/$ratio.log"
    fi
done
