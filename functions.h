#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <vector>
#include <array>
#include "vectorOperators.h"
#include <math.h>

template <typename T>
void inicialValues(T &x)
{
    for (auto& i:x)
    {
        i = (double)rand() / (4.0*RAND_MAX);
        //i=0.5;
    }
}

//double countSD(const state_type &vec);

template <typename T>
double liczWartosc(const T &vec)
{
    double result=0;
    for(auto i:vec)
    {
        result += pow(i,2);
    }
    return sqrt(result);
}

template <typename T>
void normVector(T &vec)
{
    auto val = liczWartosc(vec);
    for(auto& i:vec)
    {
        i /= val;
    }
}
//bool isSdGood(const std::vector<double>& sd, double maxVal);*/


#endif // FUNCTIONS_H_INCLUDED
