#!/bin/bash 
#$Id: fullinspect.sh 858 2007-02-12 05:46:19Z cvsrep $

#touch -t 198001010000 .early

trap '(sleep 1;[ -n "$mdirname" -a -d "$mdirname" ] && rm -vr "$mdirname"; echo killing $0; kill $$)' SIGINT

#lastdir="`ls -dt web*.dir .early |head -1`"
#potentialwebs="`ls -t web*.xml*|tac`"


#for w in `find $potentialwebs -maxdepth 0 -name 'web*.xml*' -cnewer "$lastdir"` ;do 
webnames="`find . -maxdepth 1 -name 'web*.xml*'`"
for w in `ls -1t $webnames|tac` ;do 
    echo -ne $w '\r' >&2; 
    dirname="`basename $w .gz|sed -e 's/\.bz2$//' -e 's/\.xml$/.dir/'`"
    if [ \! -d "$dirname" -o "$w" -nt "$dirname" -o \! -e "$dirname"/fixed ]; then
	webname="`pwd`/$w"
	mdirname="`pwd`/$dirname"
	mkdir -p "$dirname"
	pushd $dirname > /dev/null
	~/NewWeb/build/Inspect -a -mbipst 1e-2 -S strength.dat -Z sizes.dat -D distances.dat -B biomasses.dat -T table.dat -P img.pgm -V vul.dat -G gen.dat -L links.dat "$webname" >> inspect.log
	bzip2 --force *.dat *.log *.xml
	mdirname=""
	touch fixed
	popd > /dev/null
	touch --reference="$webname" "$dirname"
    fi
done

#rm .early

echo # leave the status line
