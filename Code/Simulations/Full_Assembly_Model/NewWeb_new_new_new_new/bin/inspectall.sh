#!/bin/bash
#$Id: inspectall.sh 627 2006-09-24 04:41:37Z cvsrep $

for w in `ls -t web*.xml*|tac`;do 
    echo -ne $w '\r' >&2; 
    stname="`basename $w .gz`"
    stname="`basename $stname .bz2`"
    stname="`basename $stname .xml`".st
    ~/NewWeb/build/Inspect -a -imspt 1e-2 -S $stname "$@" $w; 
done | \
awk '/^web_number / {n=$2} /^(T|S|Z|C|GenSD|VulSD|MxSim|Loop|Cannib|Ominv|aChnLg|aChnSD|aChnNo|aLoop|aOmniv|Ddiet|Clust|log10_bodymass|log10_plant_bodymass|log10_animal_bodymass|max_level|mean_level|level_sharpness|max_log10_bodymass|n_plants|n_animals|current_time|insertions|Zc01) / {if($2+0 == $2) print n,$2 | "sort -n > " $1 ".dat";}'

echo # leave the status line
