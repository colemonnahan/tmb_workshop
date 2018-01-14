// To demonstrate AD. A function that would be difficult to differentiate
// analytically.

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(x);
  Type y;
  y = -exp(-x*x)*pow(x+1, 2)-0.1*pow(sin(x)-3, 4.0);
  return(y);
}
