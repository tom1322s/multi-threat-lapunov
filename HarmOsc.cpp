#include "HarmOsc.h"
#include <iostream>
#include <math.h>

HarmOsc::HarmOsc()
{
    biff_type=1;
    // 1 - for k, 2 - for d, 3 - for a

    k = 2.5;
    d = 0.3;
    a = 0.1;

    updatePeriod();
}

void HarmOsc::updatePeriod()
{
    //T=2*3.14/eta;
    T = 2*3.14/sqrt(pow(d,2)-pow(a,2)); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dt=T/dtVal;
}

void HarmOsc::printInfo() // !!!!!!!!!!!!!!!!!
{
    std::cout<<std::endl;
    std::cout<<"Rozpoczynam symulacje dla\tk="<<k<<"\td="<<d<<"\ta="<<a<<std::endl;
    std::cout<<"krok : "<<dt<<std::endl;
}


bool HarmOsc::checkPoincareRequirement(const double t, const double tOld)
{
    //return POUNCARE_REQUIREMENT;
    return 1;//(sin(eta * t) > 0 && sin(eta * t)*sin(eta * tOld) < 0); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

double HarmOsc::calculateTimeToSavePoint(const double t, const double tOld)
{
    return 1;//(t - (sin(eta * t) * (t - tOld)) / (sin(eta * t) - sin(eta * tOld)));
}

void HarmOsc::changeBiffurParametr(double biffVal)
{
    switch(biff_type)
    {
    case 1:
        {
            k = biffVal;
        }
        break;

    case 2:
        {
            d = biffVal;
        }
        break;

    case 3:
        {
            a = biffVal;
        }
        break;
    }

    updatePeriod();
}

double HarmOsc::getBiffurParametr()
{
    switch(biff_type)
    {
    case 1: return k;
    case 2: return d;
    case 3: return a;

    }
    return 0;
}


void HarmOsc::operator() (const state_type &x , state_type &dxdt , const double t ) // !!!!!!!!!!!!!!!!!
{

    dxdt[0] = x[1];
    dxdt[1] = -d*x[1] -a*x[0] -pow(x[0],3) +k*(x[4]-x[0]);
    dxdt[2] = x[3];
    dxdt[3] = -d*x[3] -a*x[2] -pow(x[2],3) +k*(x[0]-x[2]);
    dxdt[4] = x[5];
    dxdt[5] = -d*x[5] -a*x[4] -pow(x[4],3) +k*(x[2]-x[4]);

}


