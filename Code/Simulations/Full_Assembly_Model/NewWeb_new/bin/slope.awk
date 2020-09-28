#!/usr/bin/awk -f
# $Id: slope.awk 802 2006-11-29 08:23:03Z cvsrep $
# Compute the slope of a list of pairs of real numbers.
# Dropps points if they are identical to predecessor, as these are
# likely to be bad data from the log_spectrum class.
{
  x=$1;
  y=$2;
  if(y!="-inf" && (NR==1 || y!=lasty)){
    sumxy+=x*y;
    sumx2+=x*x;
    sumy2+=y*y;
    sumx+=x;
    sumy+=y;
    n++;
    lasty=y;
  }
}
END{
  Sxy=(n*sumxy-sumx*sumy);
  Sxx=(n*sumx2-sumx*sumx);
  b=Sxy/Sxx;
  a=(sumy-b*sumx)/n;
  sume2=sumy2-2*b*sumxy-2*a*sumy+b*b*sumx2+2*a*b*sumx+n*a*a;
  if(n>2){
    std_b=sqrt(sume2/Sxx/(n-2));
  }else{
    std_b="inf";
  }
  print b,std_b;
}
