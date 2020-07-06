#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <vector>
#include "vectorOperators.h"

void inicialValues(state_type &x);
double countSD(const state_type &vec);
double liczWartosc(const state_type &vec);
void normVector(state_type &vec);
bool isSdGood(const std::vector<double>& sd, double maxVal);


#endif // FUNCTIONS_H_INCLUDED
