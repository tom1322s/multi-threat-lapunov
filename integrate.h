#ifndef INTEGRATE_H_INCLUDED
#define INTEGRATE_H_INCLUDED

#include <iostream>
#include <algorithm>
#include <chrono>
#include "vectorOperators.h"

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

//typedef std::vector< double > state_type;

//state_type& operator += (state_type &v1, const state_type &v2);

template < typename System,typename state_type>
void rk4(System system,const state_type &x, state_type &dx, const double t, const double dt)
{
    state_type local_x;
    std::array<state_type,4> k;
    state_type temp_k;

    local_x = x;
    system(local_x, k[0],t);
    temp_k = k[0];
    std::for_each(temp_k.begin(), temp_k.end(), [dt](double &x) {x*=0.5*dt;});
    local_x += temp_k;

    system(local_x, k[1], t+dt/2);
    temp_k = k[1];
    std::for_each(temp_k.begin(), temp_k.end(), [dt](double &x) {x*=0.5*dt;});
    local_x = x;
    local_x += temp_k;

    system(local_x, k[2],t+dt/2);
    temp_k = k[2];
    std::for_each(temp_k.begin(), temp_k.end(), [dt](double &x) {x*=dt;});
    local_x = x;
    local_x += temp_k;

    system(local_x, k[3], t+dt);

    for(unsigned int i = 0; i < x.size(); i++)
    {
        dx[i] = (k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i])*dt/6.0;
    }
}


template < typename System, typename Observer,typename state_type>
void integrateConstStare(System system,state_type &x, const double t0, const double tMax, const double dt,Observer observer)
{
    state_type dxdt;
    observer(x,t0);
    for (double t = t0; t < tMax; t+=dt)
    {
        if(t+dt<=tMax)
        {
            rk4(system,x,dxdt,t, dt);
            observer(x,t);
            x+=dxdt;
        }
        else
        {
            rk4(system,x,dxdt,t, tMax-t);
            observer(x,t);
            x+=dxdt;
        }
    }

}

template < typename System, typename Observer,typename state_type>
void integrateConst(System system,state_type &x, const double t0, const double tMax, const double dt,Observer observer)
{
    state_type dxdt;
    long int steps = (tMax-t0)/dt;
    observer(x,t0);

    for (long int i = 0; i < steps; i++)
    {
        double t = t0 + i * dt;
        rk4(system,x,dxdt,t, dt);
        observer(x,t);
        x+=dxdt;
    }


    if(double t(t0 + steps*dt);tMax - t > 0)
    {
        rk4(system,x,dxdt,t, tMax-t);
        observer(x,t);
        x+=dxdt;
    }
}


template < typename System,typename state_type>
void integrateConst(const System& system,state_type &x, const double t0, const double tMax, const double dt)
{
    state_type dxdt;
    long int steps = (tMax-t0)/dt;

    for (long int i = 0; i < steps; i++)
    {
        double t = t0 + i * dt;
        rk4(system,x,dxdt,t, dt);
        x+=dxdt;
    }


    if(double t(t0 + steps*dt);tMax - t > 0)
    {
        rk4(system,x,dxdt,t, tMax-t);
        x+=dxdt;
    }
}

template < typename System,typename state_type>
state_type integrateConstR(const System& system,state_type x, const double t0, const double tMax, const double dt)
{
    state_type dxdt;
    long int steps = (tMax-t0)/dt;

    for (long int i = 0; i < steps; i++)
    {
        double t = t0 + i * dt;
        rk4(system,x,dxdt,t, dt);
        x+=dxdt;
    }


    if(double t(t0 + steps*dt);tMax - t > 0)
    {
        rk4(system,x,dxdt,t, tMax-t);
        x+=dxdt;
    }
    return x;
}

template < typename System,typename state_type>
state_type integrateConstRStare(const System& system,state_type x, const double t0, const double tMax, const double dt)
{
    state_type dxdt;

    for (double t = t0; t < tMax; t+=dt)
    {
        if(t+dt<=tMax)
        {
            rk4(system,x,dxdt,t, dt);
            x+=dxdt;
        }
        else
        {
            rk4(system,x,dxdt,t, tMax-t);
            x+=dxdt;
        }
    }

    return x;
}

/* !!!!template < typename System>
state_type integrateConst(const System& system,state_type x, const double t0, const double tMax, const double dt)
{
const auto startTime = std::chrono::high_resolution_clock::now();
    state_type dxdt(x.size());
    std::cout<<tMax-t0<<std::endl;
    int counter = 0;
    // zmienic t z double na ilosc ktokow ile chce wykonac liczyc kolejne czas jako dt razy licznik pêtli
    for (double t = t0; t < tMax; t+=dt)
    {
        counter++;
        //if(t+dt<=tMax)
        {

            rk4(system,x,dxdt,t, dt);

            x+=dxdt;

        }
        //else
        //{
          //  rk4(system,x,dxdt,t, tMax-t);
            //x+=dxdt;
        //}
    }
    std::cout<<counter<<std::endl;
    const auto endTime = std::chrono::high_resolution_clock::now();

    std::cout<<(duration_cast<duration<double, std::milli>>(endTime-startTime).count())<<std::endl;

    return x;
}*/



#endif // INTEGRATE_H_INCLUDED
