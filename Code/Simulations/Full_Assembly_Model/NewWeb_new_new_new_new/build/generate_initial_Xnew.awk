#!/usr/bin/awk -f
# $Id: generate_initial_Xnew.awk 234 2006-02-02 06:03:12Z cvsrep $
# extracts the lines containing data in space matrix file

hrules==2 {
  line++;
  Variable[line]=$2
  Min[line]=$3
  Max[line]=$4
}
/^--/ {hrules++;}
END{
  srand();
  of = "Xnew.mat";
  print "X matrix" > of;
  print "------------------------------------" > of;
  printf "Case " > of;
  for(i=1;i<=line;i++){
    printf "%s ", Variable[i] > of;
  }
  print "" > of;
  print "------------------------------------" > of;
  for(j=1;j<=starting_lines;j++){
    printf "%i ", j > of;
    for(i=1;i<=line;i++){
      printf "%s ", Min[i]+rand()*(Max[i]-Min[i]) > of;
    }
    print "" > of;
  }
} 
