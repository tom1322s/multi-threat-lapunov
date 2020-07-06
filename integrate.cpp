#include "integrate.h"

state_type& operator += (state_type &v1, const state_type &v2)
{

    for(unsigned int i=0;i<v1.size();i++)
    {
        v1[i]+=v2[i];
    }
    return v1;
}
