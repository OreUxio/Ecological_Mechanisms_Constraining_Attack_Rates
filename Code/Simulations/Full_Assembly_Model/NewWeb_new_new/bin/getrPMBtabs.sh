#!/bin/bash
#$Id: getrPMBtabs.sh 795 2006-11-28 23:55:24Z cvsrep $
for w in `ls -t web*.xml*|tac`;do 
    echo -ne $w '\r' >&2; 
    MBname="`basename $w .gz`"
    MBname="`basename $MBname .bz2`"
    MBname="`basename $MBname .xml`"
    ~/NewWeb/build/NewWorld $w ~/NewWeb/build/print_each_step.cfg > /dev/null
    ~/NewWeb/build/NewWorld relaxed.xml ~/NewWeb/build/print_each_step.cfg > /dev/null
    ~/NewWeb/build/Inspect -a -imst 1e-2 -T tmp.tab relaxed.xml|grep n_plants ;
    awk 'int($5+0.5)==1 {print $2,$3,int($5+1.5);}' tmp.tab | bzip2 > $MBname.PMB.bz2
done 
rm tmp.tab

echo >&2 # leave the status line
