#!/bin/bash

for f in *; do
    if [ -d "$f" ]; then
    # $f is a directory
       ratio=${f#r_*}
       echo $ratio
       ../bin/fastVQA --vqe -e -r 15 -x $ratio &> "log/$ratio.log"
    fi
done
