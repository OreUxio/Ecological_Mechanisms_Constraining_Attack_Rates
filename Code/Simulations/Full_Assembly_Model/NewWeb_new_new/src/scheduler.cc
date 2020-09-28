// $Id: scheduler.cc 371 2006-03-31 22:40:51Z cvsrep $

#include "scheduler.h"

const double repeat_after_t::not_set=-1;

double scheduler_t::
time_of_next_call() const{
  double t_min=HUGE_VAL;
  for(scheduler_entry_lt::const_iterator i = the_entries.begin();
      i!=the_entries.end();
      i++){
    if((*i)->is_active()){
      double nc=(*i)->next_call();
      if(nc<t_min)
	t_min=nc;
    }
  }
  return t_min;
}


scheduler_t::scheduler_t():the_entries(){};
scheduler_t & scheduler_t::add(scheduler_entry_t & e){
  the_entries.push_back(&e);
  return *this;
}
scheduler_t & scheduler_t::remove(scheduler_entry_t & e){
  the_entries.remove(&e);
  return *this;
}

scheduler_entry_t::~scheduler_entry_t(){};
bool scheduler_entry_t::is_active(){
  return this_is_active;
}

double repeat_after_t::next_call() const {
  return the_next_call;
}

repeat_after_t::repeat_after_t(double dt,double start_time) : 
  the_dt(dt?dt:HUGE_VAL), 
  the_next_call(start_time), 
  the_time_since_last_call(not_set)
{
  this_is_active=(dt>0);
};

bool repeat_after_t::due_now(double t) {
  if(this_is_active && t>=next_call() ){
    the_time_since_last_call=t-the_next_call+the_dt;
    the_next_call=t+the_dt;
    return true;
  }else{
    return false;
  }
}

double repeat_after_t::time_since_last_call() const{
  if(the_time_since_last_call!=not_set)
    return the_time_since_last_call;
  else
    FATAL_ERROR("time since last call cannot be determined");
  return 0;
}

void repeat_after_t::reset(double t){
  the_next_call=t;
  the_time_since_last_call=not_set;
}

do_after_every_t::do_after_every_t(double dt) :
  repeat_after_t(dt,dt) {
}

at_multiples_of_t::at_multiples_of_t(double dt,double start_time) :
  repeat_after_t(dt,((long int)(start_time/dt)+1)*dt) {
}

bool at_multiples_of_t::due_now(double t){
  bool hold=repeat_after_t::due_now(t);
  if(hold){
    //correct rounding errors:
    the_next_call=(long int)(the_next_call/the_dt+0.5)*the_dt;
  }
  return hold;
}

void at_multiples_of_t::reset(double t){
  repeat_after_t::reset((long int)(t/the_dt)+1);
}
