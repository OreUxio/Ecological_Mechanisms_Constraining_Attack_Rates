#!/usr/bin/awk -f
# $Id: monotones.awk 354 2006-03-23 14:28:56Z cvsrep $
# extract the size of variations over monotone ups and down in time series

BEGIN {
  going_up=1;
  going_down=2;
  unclear=3;
  going_one_up=4;
  state=unclear;
}

{
#(debug): print state, $1, last ,last_start;
 if(state==unclear){
    last_start=$1;
    state=going_up; #guessing
  }else if(state==going_up){
    if($1<last){
#passed maximum
      print last-last_start;
      last_start=last;
      state=going_down;
    }
  }else if(state==going_down){
    if($1>last){
#passed minimum
      if($1-last==1){
	state=going_one_up;
      }else{
	print last-last_start;
	last_start=last;
	state=going_up;
      }
    }
  }else if(state==going_one_up){
    if($1==last){
# keep waiting what happens
      state=going_one_up;
    }else if($1>last){
# we did not fall back, extinction wave has stopped.
      print last-1-last_start;
      last_start=last-1;
      state=going_up;
    }else{
      state=going_down;
    }
  }
  else {
    print "unknown state!!" > /dev/stderr;
    exit(1);
  }
  last=$1;
}

