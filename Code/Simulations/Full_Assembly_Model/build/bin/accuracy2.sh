#!/bin/bash

for f in fitness?*.dat; do
    step_count="`head -1 $f|awk '{print $1}'`"
    webname="`echo $f | sed 's/fitness/CompetiWeb/'`"
    echo -n "$step_count " >&2
    acc=`~/NewWeb/build/CompetiWeb accuracy2 "$webname" | awk '/a1o.rms_max_diff...mean/ {print $3}'`
    [ $? == 0 ] || (echo >&2 ; echo Problem with $webname >&2)
    echo $step_count $acc
done |\
awk '{if($2 != "")print;}'| sort -n > accuracy2.dat

echo >&2