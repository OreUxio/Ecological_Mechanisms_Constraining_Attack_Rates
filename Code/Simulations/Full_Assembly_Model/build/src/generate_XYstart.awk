#!/usr/bin/awk -f
# $Id: generate_XYstart.awk 234 2006-02-02 06:03:12Z cvsrep $
# extracts the lines containing data in space matrix file

hrules==2 {
  line++;
  Variable[line]=$2
  Min[line]=$3
  Max[line]=$4
}
/^--/ {hrules++;}
END{
  of = "Xstart.mat";
  print "X matrix" > of;
  print "------------------------------------" > of;
  printf "Case " > of;
  for(i=1;i<=line;i++){
    printf "%s ", Variable[i] > of;
  }
  print "" > of;
  print "------------------------------------" > of;
  
  of = "Ystart.mat";
  print "------------------------------" > of;
  print "Case y" > of;
  print "------------------------------" > of;
} 
