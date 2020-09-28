#!/bin/csh
#$Id: b1.bat 2502 2017-02-27 17:35:59Z axel $

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
