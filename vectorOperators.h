#ifndef VECTOROPERATORS_H_INCLUDED
#define VECTOROPERATORS_H_INCLUDED

#include <stdlib.h>
#include <array>


template <typename state_type>
double operator & (const state_type& lhs, const state_type& rhs)
{
    double result = 0;
    for (size_t i = 0; i < lhs.size(); i++)
    {
        result += lhs[i] * rhs[i];
    }
    return result;
}

template <typename state_type>
state_type operator * (const state_type& lhs, const state_type& rhs)
{
    state_type result;
    for (size_t i = 0; i < lhs.size(); i++)
    {
        result[i] = lhs[i] * rhs[i];
    }
    return result;
}

template <typename state_type>
state_type odj (const state_type& lhs, const state_type& rhs)
{
    state_type result;
    for (size_t i = 0; i < lhs.size(); i++)
    {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

template <typename state_type>
state_type operator * (const state_type& scr, double scalar)
{
    state_type result;
    for(unsigned int i=0; i < scr.size(); i++)
    {
        result[i]=(scr[i] * scalar);
    }
    return result;
}


template <typename state_type>
state_type& operator += (state_type &v1, const state_type &v2)
{

    for(unsigned int i=0;i<v1.size();i++)
    {
        v1[i]+=v2[i];
    }
    return v1;
}

template <typename state_type>
state_type& operator -= (state_type &vec, const state_type &scr)
{
    if(vec.size() != scr.size()) abort();
    for(unsigned int i=0; i < vec.size(); i++)
    {
        vec[i] -= scr[i];
    }
    return vec;
}

template <typename state_type>
state_type& operator /= (state_type &vec, double x)
{
    for(auto& i:vec)
    {
        i /= x;
    }
    return vec;
}
#endif // VECTOROPERATORS_H_INCLUDED
