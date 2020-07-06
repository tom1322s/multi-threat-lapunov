#ifndef INTEGRATE_H_INCLUDED
#define INTEGRATE_H_INCLUDED

#include <iostream>
#include <algorithm>
#include <chrono>

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

typedef std::vector< double > state_type;

state_type& operator += (state_type &v1, const state_type &v2);

template < typename System>
void rk4(System system,const state_type &x, state_type &dx, const double t, const double dt)
{
    state_type local_x;
    state_type k0(x.size());
    state_type k1(x.size());
    state_type k2(x.size());
    state_type k3(x.size());
    state_type temp_k;

    local_x = x;
    system(local_x, k0,t);
    temp_k = k0;
    for_each(temp_k.begin(), temp_k.end(), [dt](double &x) {x*=0.5*dt;});
    local_x += temp_k;

    system(local_x, k1, t+dt/2);
    temp_k = k1;
    for_each(temp_k.begin(), temp_k.end(), [dt](double &x) {x*=0.5*dt;});
    local_x = x;
    local_x += temp_k;

    system(local_x, k2,t+dt/2);
    temp_k = k2;
    for_each(temp_k.begin(), temp_k.end(), [dt](double &x) {x*=dt;});
    local_x = x;
    local_x += temp_k;

    system(local_x, k3, t+dt);

    for(unsigned int i = 0; i < x.size(); i++)
    {
        dx[i] = (k0[i] + 2*k1[i] + 2*k2[i] + k3[i])*dt/6.0;
    }
}

template < typename System, typename Observer>
void integrateConst(System system,state_type &x, const double t0, const double tMax, const double dt,Observer observer)
{
    state_type dxdt(x.size());
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

/*template < typename System>
void integrateConst(const System& system,state_type &x, const double t0, const double tMax, const double dt)
{
    state_type dxdt(x.size());
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

}*/

template < typename System>
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
        /*else
        {
            rk4(system,x,dxdt,t, tMax-t);
            x+=dxdt;
        }*/
    }
    std::cout<<counter<<std::endl;
    const auto endTime = std::chrono::high_resolution_clock::now();

    std::cout<<(duration_cast<duration<double, std::milli>>(endTime-startTime).count())<<std::endl;

    return x;
}



#endif // INTEGRATE_H_INCLUDED
