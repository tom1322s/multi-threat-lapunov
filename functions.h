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
template <typename T>
T countSD(const std::vector<T> &vec)
{
    T averege;
    for(auto& i:averege) i=0.0;
    for_each(vec.begin(),vec.end(),[&averege](T i) {averege+=i;});
    averege/=vec.size();

    T sd;
    for(auto& i:sd) i=0.0;
    for_each(vec.begin(),vec.end(),[&sd,averege](T i) {sd+=odj(i,averege)*odj(i,averege);});
    sd/=vec.size();
    return sd;
}

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

template <typename T>
bool isSdGood(const T& sd, double maxVal)
{
    bool result = false;
    for(auto i:sd)
    {
        result |= (i > maxVal);
    }
    return result;
}


#endif // FUNCTIONS_H_INCLUDED
