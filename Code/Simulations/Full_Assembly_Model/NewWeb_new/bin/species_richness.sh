#!/bin/bash

for f in fitness?*.dat; do
    step_count="`head -1 $f|awk '{print $1}'`"
    webname="`echo $f | sed 's/fitness/CompetiWeb/'`"
    S=`~/NewWeb/build/CompetiWeb eval "$webname" | awk '/number_of_spe/ {print $3}'`
    [ $? == 0 ] || echo $webname >&2
    echo $step_count $S
done |\
awk '{if($2 != "")print;}'| sort -n > richness.dat
