#include "vectorOperators.h"

double operator & (const state_type& lhs, const state_type& rhs)
{
    double result = 0;
    for (int i = 0; i < lhs.size(); i++)
    {
        result += lhs[i] * rhs[i];
    }
    return result;
}

state_type operator * (const state_type& scr, double scalar)
{
    state_type result;
    result.reserve(scr.size());
    for(unsigned int i=0; i < scr.size(); i++)
    {
        result.push_back(scr[i] * scalar);
    }
    return result;
}


state_type& operator -= (state_type &vec, const state_type &scr)
{
    if(vec.size() != scr.size()) abort();
    for(unsigned int i=0; i < vec.size(); i++)
    {
        vec[i] -= scr[i];
    }
    return vec;
}
