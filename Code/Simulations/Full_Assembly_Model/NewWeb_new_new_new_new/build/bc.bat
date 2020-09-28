#!/bin/csh
#$Id: bc.bat 2502 2017-02-27 17:35:59Z axel $

cd /hwork0/d6/ishii/${QUEUENAME}

setenv LANG C
setenv C_PROGINF DETAIL
setenv F_FTRACE YES
echo $SHELL
printenv

hostname
date

cat > run.cfg <<EOF
observation_time=(10^10)*years

EOF

./NewWorld SAVED.xml NewWorld.cfg > log.dat

date

cat log.dat >> NewWorld.dat
cat SAVED.xml >> SAVED.xml.bac

echo batch job finished
