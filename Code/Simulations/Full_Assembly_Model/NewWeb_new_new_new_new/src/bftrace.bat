#!/bin/csh
#$Id: bftrace.bat 2502 2017-02-27 17:35:59Z axel $

cd /hwork0/d6/ishii/${QUEUENAME}

setenv LANG C
setenv C_PROGINF DETAIL
setenv F_FTRACE YES
echo $SHELL
printenv

hostname
date

setenv currentTime 111751.44073182516 #`awk -F '[>]\|[<]' '/current_time/{print $3;exit}' START.xml`

cat > run.cfg <<EOF
observation_time=(${currentTime} + 60.1)*years
time_between_printout=1*year
time_between_analysis=0
TRACEFLAG=(1+2+4+8+16+32)
EOF

time ./NewWorld START.xml run.cfg > NewWorld.dat

date

sxftrace++ > ftrace.txt

date

echo batch job finished
