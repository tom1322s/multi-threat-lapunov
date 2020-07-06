#include "functions.h"
#include <algorithm>



/*double countSD(const state_type &vec)
{
    double averege=0;
    for_each(vec.begin(),vec.end(),[&averege](double i) {averege+=i;});
    averege/=vec.size();

    double sd=0;
    for_each(vec.begin(),vec.end(),[&sd,averege](double i) {sd+=(i-averege)*(i-averege);});
    sd/=vec.size();
    return sd;
}

double liczWartosc(const state_type &vec)
{
    double result=0;
    for(auto i:vec)
    {
        result += pow(i,2);
    }
    return sqrt(result);
}

void normVector(state_type &vec)
{
    auto val = liczWartosc(vec);
    for(auto& i:vec)
    {
        i /= val;
    }
}

bool isSdGood(const std::vector<double>& sd, double maxVal)
{
    bool result = false;
    for(auto i:sd)
    {
        result |= (i > maxVal);
    }
    return result;
}*/
