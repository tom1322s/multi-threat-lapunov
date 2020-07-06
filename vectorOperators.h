#ifndef VECTOROPERATORS_H_INCLUDED
#define VECTOROPERATORS_H_INCLUDED

#include <stdlib.h>
#include <vector>

typedef std::vector< double > state_type;


double operator & (const state_type& lhs, const state_type& rhs);
state_type operator * (const state_type& scr, double scalar);
state_type& operator -= (state_type &vec, const state_type &scr);
#endif // VECTOROPERATORS_H_INCLUDED
