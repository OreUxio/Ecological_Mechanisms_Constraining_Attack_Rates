#!/bin/bash
#$Id: getrMBtab.sh 807 2006-11-29 08:36:28Z cvsrep $
for w in `ls -t web*.xml*|tac`;do 
    echo -ne $w '\r' >&2; 
    MBname="`basename $w .gz`"
    MBname="`basename $MBname .bz2`"
    MBname="`basename $MBname .xml`"
    ~/NewWeb/build/NewWorld $w ~/NewWeb/build/print_each_step.cfg > /dev/null 2>&1
    ~/NewWeb/build/NewWorld relaxed.xml ~/NewWeb/build/print_each_step.cfg > /dev/null 2>&1
    ~/NewWeb/build/NewWorld relaxed.xml ~/NewWeb/build/print_each_step.cfg > /dev/null 2>&1
    ~/NewWeb/build/NewWorld relaxed.xml ~/NewWeb/build/print_each_step.cfg > /dev/null 2>&1
    ~/NewWeb/build/Inspect -a -it 1e-2 -T tmp.tab relaxed.xml>/dev/null ;
    awk '{print $2,$3,int($5+1.5);}' tmp.tab | bzip2 > $MBname.MB.bz2
done 
rm tmp.tab

echo >&2 # leave the status line
