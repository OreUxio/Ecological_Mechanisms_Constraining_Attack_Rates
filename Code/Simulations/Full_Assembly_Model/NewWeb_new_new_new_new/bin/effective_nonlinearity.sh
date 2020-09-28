#!/bin/bash

for f in fitness?*.dat; do
    step_count="`head -1 $f|awk '{print $1}'`"
    webname="`echo $f | sed 's/fitness/CompetiWeb/'`"
    ~/NewWeb/build/CompetiWeb logistic "$webname" >&2 #> /dev/null 2>&1
    awk -v c=$step_count '{print c,$2/$1}' logistic.dat
done |\
sort -n > nonlinearity.dat
