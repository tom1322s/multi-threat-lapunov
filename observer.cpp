#include "observer.h"
/*
double liczWartosc(state_type &x, int numberPer)
{
    double suma = 0;

    for (int i = ORDER * numberPer; i < ORDER * (numberPer + 1); i++)
    {
        suma += pow(x[i],2);
    }
    return sqrt(suma);
}

void normVector(state_type &x, int numberPer)
{
    double module = liczWartosc(x, numberPer);
    for (int i = ORDER * numberPer; i < ORDER * (numberPer + 1); i++)
    {
        x[i] /= module;
    }
}

void normPerturbationAndUpdatePerpendicularComponent(state_type &x)
{
    state_type z1(ORDER), z2(ORDER);
    copyPerturbation(z1, z2, x);

    z2 -= z1 * ((z1 & z2) / (z1 & z1));


    for(unsigned int i = ORDER * 2; i < O2ORDER; i++)
    {
        x[i] = z2[i - ORDER * 2];
    }

    double n1, n2;
    n1 = liczWartosc(x,1);
    n2 = liczWartosc(x,2);
    if (n1 < pow(10,-6) || n1 > pow(10,6)) normVector(x,1);
    if (n2 < pow(10,-6) || n2 > pow(10,6)) normVector(x,2);
}

void copyPerturbation(state_type &vec1, state_type &vec2, const state_type &scr)
{
    for(unsigned int i = ORDER; i <ORDER * 2; i++)
    {
        vec1[i - ORDER] = scr[i];
    }
    for(unsigned int i = ORDER * 2; i <O2ORDER; i++)
    {
        vec2[i - ORDER * 2] = scr[i];
    }
}

double operator & (const state_type &vec1, const state_type &vec2)
{
    double result = 0;
    if(vec1.size() != vec2.size()) abort();

    for(unsigned int i = 0; i < vec1.size(); i++)
    {
        result += vec1[i] * vec2[i];
    }
    return result;
}

state_type operator * (state_type& scr, double scalar)
{
    state_type result;
    result.reserve(ORDER);
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
*/


