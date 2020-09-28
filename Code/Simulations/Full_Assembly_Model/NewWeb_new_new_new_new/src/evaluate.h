/* $Id: evaluate.h 2250 2013-01-04 20:16:04Z axel $ */

#ifndef __EVALUATE_H__
#define __EVALUATE_H__

#ifdef __cplusplus
extern "C" {
#endif

  int evaluate_expression(char ** value);
  void set_evaluator_variable(const char * name,double value);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* define some origian C++ stuff, too: */
#include "CLHEP_Evaluator.h"

int CXXevaluate_expression(char ** value);

/// Encapsulates a formula evaluator.
class my_evaluator_t : private HepTool::Evaluator {
  int f;
  bool if_error_hint_at_position(const char * const & value);
 public:
  my_evaluator_t();
  double operator()(const char * const & value);
  void set_variable(const char * name, double value);

  friend bool if_error_hint_at_position(char * const &  value,my_evaluator_t & eval);
  friend int evaluate_expression(char ** value);
  friend void set_evaluator_variable(const char * name, double value);
};

extern my_evaluator_t eval;    

#endif

#endif /* __EVALUATE_H__ */
