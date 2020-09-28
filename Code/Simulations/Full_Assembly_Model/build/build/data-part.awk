#!/usr/bin/awk -f
# $Id: data-part.awk 234 2006-02-02 06:03:12Z cvsrep $
# extracts the lines containing data in space matrix file

hrules==2 {print;}
/^--/ {hrules++;}
