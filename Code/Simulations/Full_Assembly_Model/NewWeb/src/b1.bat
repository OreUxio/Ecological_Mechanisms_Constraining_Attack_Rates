#!/bin/csh
#$Id: b1.bat 376 2006-04-01 12:30:37Z cvsrep $

cd /hwork0/d6/ishii/${QUEUENAME}

setenv LANG C
setenv C_PROGINF DETAIL
echo $SHELL
printenv

hostname
date

time ./NewWorld NewWorld.cfg > NewWorld.dat

date

echo batch job finished
