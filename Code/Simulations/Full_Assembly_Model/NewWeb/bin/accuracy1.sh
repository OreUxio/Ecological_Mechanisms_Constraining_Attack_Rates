#!/bin/bash

for f in fitness?*.dat; do
    step_count="`head -1 $f|awk '{print $1}'`"
    webname="`echo $f | sed 's/fitness/CompetiWeb/'`"
    echo -n "$step_count " >&2
    acc="`~/NewWeb/build/CompetiWeb accuracy1 "$webname" | awk '/a1o.rms_max_diff...mean/ {raw=$3} /a1o.rms_max_diff_corrected...mean/ {corrected=$3} END {print raw, corrected}'`"
    [ $? == 0 ] || (echo >&2 ; echo Problem with $webname >&2)
    echo $step_count $acc
done |\
awk '{if($2 != "")print;}'| sort -n > accuracy1.dat

echo >&2