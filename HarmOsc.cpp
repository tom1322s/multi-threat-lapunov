#include "HarmOsc.h"
#include <iostream>
#include <math.h>


#if OSC_KIND == 0
HarmOsc::HarmOsc()
{
    biff_type=1;
    // 1 - for k, 2 - for d, 3 - for a

    k = 0;
    d = 0.3;
    a = 0.1;

    updatePeriod();
}

void HarmOsc::updatePeriod()
{
    //T=2*3.14/eta;
    //T = 2*3.14/sqrt(pow(d,2)-pow(a,2)); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //dt=T/dtVal;
    dt = 0.0005;//0.0002;//0.0005;
    T =1.0;

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
//2
    //dxdt[0] = x[1];
    //dxdt[1] = -d*x[1] -a*x[0] -pow(x[0],3);


    dxdt[0] = x[1];
    dxdt[1] = -d*x[1] -a*x[0] -pow(x[0],3) +k*(x[4]-x[0]);
    dxdt[2] = x[3];
    dxdt[3] = -d*x[3] -a*x[2] -pow(x[2],3) +k*(x[0]-x[2]);
    dxdt[4] = x[5];
    dxdt[5] = -d*x[5] -a*x[4] -pow(x[4],3) +k*(x[2]-x[4]);

    /*dxdt[0] = x[1];
    dxdt[1] = -d*x[1] -a*x[0] -pow(x[0],3) +k*(x[16]-x[0]);
    dxdt[2] = x[3];
    dxdt[3] = -d*x[3] -a*x[2] -pow(x[2],3) +k*(x[0]-x[2]);
    dxdt[4] = x[5];
    dxdt[5] = -d*x[5] -a*x[4] -pow(x[4],3) +k*(x[2]-x[4]);
    dxdt[6] = x[7];
    dxdt[7] = -d*x[7] -a*x[6] -pow(x[6],3) +k*(x[4]-x[6]);
    dxdt[8] = x[9];
    dxdt[9] = -d*x[9] -a*x[8] -pow(x[8],3) +k*(x[6]-x[8]);
    dxdt[10] = x[11];
    dxdt[11] = -d*x[11] -a*x[10] -pow(x[10],3) +k*(x[8]-x[10]);
    dxdt[12] = x[13];
    dxdt[13] = -d*x[13] -a*x[12] -pow(x[12],3) +k*(x[10]-x[12]);
    dxdt[14] = x[15];
    dxdt[15] = -d*x[15] -a*x[14] -pow(x[14],3) +k*(x[12]-x[14]);
    dxdt[16] = x[17];
    dxdt[17] = -d*x[17] -a*x[16] -pow(x[16],3) +k*(x[14]-x[16]);
    dxdt[18] = x[19];
    dxdt[19] = -d*x[19] -a*x[18] -pow(x[18],3) +k*(x[16]-x[18]);
    dxdt[20] = x[21];
    dxdt[21] = -d*x[21] -a*x[20] -pow(x[20],3) +k*(x[18]-x[20]);*/

}
#elif OSC_KIND == 1
HarmOsc::HarmOsc()
{
    biff_type=2;
    // 1 - for h, 2 - for q, 3 - for eta

    h = 0.05;
    q = 1.1;
    eta = 0.33;

    updatePeriod();
}

void HarmOsc::updatePeriod()
{
    //T=2*3.14/eta;
    //T = 2*3.14/sqrt(pow(d,2)-pow(a,2)); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //dt=T/dtVal;
    dt = 0.0005;//0.0008;
    T =6.0;

}

void HarmOsc::printInfo() // !!!!!!!!!!!!!!!!!
{
    std::cout<<std::endl;
    std::cout<<"Rozpoczynam symulacje dla\th="<<h<<"\tq="<<q<<"\teta="<<eta<<std::endl;
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
            h = biffVal;
        }
        break;

    case 2:
        {
            q = biffVal;
        }
        break;

    case 3:
        {
            eta = biffVal;
        }
        break;
    }

    updatePeriod();
}

double HarmOsc::getBiffurParametr()
{
    switch(biff_type)
    {
    case 1: return h;
    case 2: return q;
    case 3: return eta;

    }
    return 0;
}


void HarmOsc::operator() (const state_type &x , state_type &dxdt , const double t ) // !!!!!!!!!!!!!!!!!
{

    dxdt[0] = x[1];
    dxdt[1] = -h*x[1] -pow(x[0],3) +q*sin(eta*t);
    dxdt[2] = x[3];
    dxdt[3] = -h*x[3] -pow(x[2],3) +q*sin(eta*t);
    dxdt[4] = x[5];
    dxdt[5] = -h*x[5] -pow(x[4],3) +q*sin(eta*t);
    /*dxdt[6] = x[7];
    dxdt[7] = -h*x[7] -pow(x[6],3) +q*sin(eta*t);
    dxdt[8] = x[9];
    dxdt[9] = -h*x[9] -pow(x[8],3) +q*sin(eta*t);
    dxdt[10] = x[11];
    dxdt[11] = -h*x[11] -pow(x[10],3) +q*sin(eta*t);
    dxdt[12] = x[13];
    dxdt[13] = -h*x[13] -pow(x[12],3) +q*sin(eta*t);
    dxdt[14] = x[15];
    dxdt[15] = -h*x[15] -pow(x[14],3) +q*sin(eta*t);
    dxdt[16] = x[17];
    dxdt[17] = -h*x[17] -pow(x[16],3) +q*sin(eta*t);
    dxdt[18] = x[19];
    dxdt[19] = -h*x[19] -pow(x[18],3) +q*sin(eta*t);
    dxdt[20] = x[21];
    dxdt[21] = -h*x[21] -pow(x[20],3) +q*sin(eta*t);*/

}
#elif OSC_KIND == 2

HarmOsc::HarmOsc()
{
    biff_type=3;
    // 1 - for ksiY, 2 - for q, 3 - for eta, 4 - for ni, 5 - for detla, 6 - for ksi

    ksiY= 0.0884;
    q= 0.0769;
    eta= 0.95;
    ni= 0.6250/((ORDER-2)/2);
    delta= 0.0604;
    ksi= 0.1116;

    updatePeriod();
}

void HarmOsc::updatePeriod()
{
    //T=2*3.14/eta;
    //T = 2*3.14/sqrt(pow(d,2)-pow(a,2)); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //dt=T/dtVal;
    dt = 0.0005;
    T =6.68;

}

void HarmOsc::printInfo() // !!!!!!!!!!!!!!!!!
{
    std::cout<<std::endl;
    std::cout<<"Rozpoczynam symulacje dla\tksiY="<<ksiY<<"\tq="<<q<<"\teta="<<eta<<"\tni="<<ni<<"\tdelta="<<delta<<"\tksi="<<ksi<<std::endl;
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
            ksiY = biffVal;
        }
        break;

    case 2:
        {
            q = biffVal;
        }
        break;

    case 3:
        {
            eta = biffVal;
        }
        break;
    case 4:
        {
            ni = biffVal;
        }
        break;

    case 5:
        {
            delta = biffVal;
        }
        break;

    case 6:
        {
            ksi = biffVal;
        }
        break;
    }

    updatePeriod();
}

double HarmOsc::getBiffurParametr()
{
    switch(biff_type)
    {
    case 1: return ksiY;
    case 2: return q;
    case 3: return eta;
    case 4: return ni;
    case 5: return delta;
    case 6: return ksi;

    }
    return 0;
}


void HarmOsc::operator() (const state_type &x , state_type &dxdt , const double t ) // !!!!!!!!!!!!!!!!!
{
//2
    //dxdt[0] = x[1];
    //dxdt[1] = -ksiY*x[1] -x[0] +q*sin(eta*t);

//4+
    double a=0,b=0;
    for(int i = 1; i<ORDER/2;i++)
    {
        a += pow(sin(x[i*2]),2)*delta*ni+sin(x[i*2])*x[i*2+1]*ksi*ni-cos(x[i*2])*pow(x[i*2+1],2)*ni;
        b += ni*sin(x[i*2])*sin(x[i*2]);
    }



    dxdt[0] = x[1];
    dxdt[1] = -(a+q*sin(eta*t)-ksiY*x[1]-x[0])/(b-1);

    if(b-1<10e-6 && b-1>-10e-6)
    {
        std::cout<<"blad"<<std::endl;
        dxdt[1] = 0;
    }

    for(int i = 1; i<ORDER/2;i++)
    {
        dxdt[i*2] = x[i*2+1];
        dxdt[i*2+1] = -(dxdt[1]+delta)*sin(x[i*2]) -ksi*x[i*2+1];
    }





}

#elif OSC_KIND == 3

HarmOsc::HarmOsc()
{
    biff_type=3;
    // 1 - for alfa, 2 - for beta, 3 - for k

    alfa = 0.1;
    beta = 1;
    k = 1;

    updatePeriod();
}

void HarmOsc::updatePeriod()
{
    //T=2*3.14/eta;
    //T = 2*3.14/sqrt(pow(d,2)-pow(a,2)); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //dt=T/dtVal;
    //T =1.4;
    T = 2*3.14/4.6403;
    dt =T/37890.0;

}

void HarmOsc::printInfo() // !!!!!!!!!!!!!!!!!
{
    std::cout<<std::endl;
    std::cout<<"Rozpoczynam symulacje dla\talfa="<<alfa<<"\tbeta="<<beta<<"\tk="<<k<<std::endl;
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
            alfa = biffVal;
        }
        break;

    case 2:
        {
            beta = biffVal;
        }
        break;

    case 3:
        {
            k = biffVal;
        }
        break;

    }

    updatePeriod();
}

double HarmOsc::getBiffurParametr()
{
    switch(biff_type)
    {
    case 1: return alfa;
    case 2: return beta;
    case 3: return k;

    }
    return 0;
}


void HarmOsc::operator() (const state_type &x , state_type &dxdt , const double t ) // !!!!!!!!!!!!!!!!!
{
//2
    //dxdt[0] = x[1];
    //dxdt[1] = alfa*(1-x[0]*x[0])*x[1]-beta*x[0];

//4+

    for(int i = 0; i<ORDER;i+=2)
    {
        double xPrev;
        if(i==0) xPrev = x[ORDER-2];
        else xPrev = x[i-2];

        dxdt[i] = x[i+1];
        dxdt[i+1] = alfa*(1-x[i]*x[i])*x[i+1]-beta*x[i] + k* (xPrev - x[i]);// + 12.95*sin(4.6403*t);
    }




}

#endif // OSC_KIND

