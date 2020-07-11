#include "programs.h"
#include "functions.h"
#include <thread>
#include <future>
#include <fstream>
//#include <math.h>
#include <chrono>

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;




void timePhraseLaypunov()
{
    std::fstream FileL, FileT, FileP;
    FileP.open( "C:\\chwilowy\\prog\\timePhrase.txt", std::ios::out );
    FileL.open( "C:\\chwilowy\\prog\\timeLaponov.txt", std::ios::out );
    FileT.open( "C:\\chwilowy\\prog\\time.txt", std::ios::out );
    if( FileP.good() == false ||  FileT.good() == false || FileL.good() == false)
    {
        std::cout<<"blad"<<std::endl;
        //std::cerr << "Error: " << strerror(errno);
        exit(250);
    }

    HarmOsc Osc;
    std::vector<double> times;
    std::vector<std::vector<double>> x_vec;
    std::vector<std::vector<double>> lambda;
    pushBackTimeAndState<HarmOsc> bufferInfo( times, x_vec,Osc);

    state_type x;
    inicialValues(x);

    std::array<state_type,ORDER> spectrum_z;
    for(auto& i:spectrum_z)
    {
        inicialValues(i);
    }

    Osc.printInfo();
    double delta = 1e-6;
    state_type lambdaSum;
    std::vector<std::vector<double>> calTime;

    double t=0.0;

#if TIME_SD
    state_type sd;
    for(auto& i:sd) i = 1.0;
    std::vector<state_type> lambdaBuffer;
#endif // TIME_SD

#if TIME_SD
    for (long int step = 0; step < Osc.nouberOfPeriodSkiped || isSdGood(sd, Osc.sdVal); step++)
#else
    for (long int step = 0; step < Osc.nouberOfPeriod; step++)
#endif // TIME_SD
    {
        t = Osc.t0 + step * Osc.T;
        const auto startTime = std::chrono::high_resolution_clock::now();


        std::array<state_type,ORDER> v;

#if MULTI_THREAD

        std::array<std::future<state_type>,ORDER> thread;
        for (int i = 0; i < Osc.order; i++)
        {
            v[i]=x;
            v[i][i]+=delta;
            thread[i] = std::async(std::launch::async, &integrateConstR<HarmOsc,state_type>, Osc,v[i],t,t+Osc.T,Osc.dt);
        }
        integrateConst(Osc,x,t,t+Osc.T,Osc.dt,bufferInfo);
        for (int i = 0; i < Osc.order; i++)
        {
            v[i] = thread[i].get();
        }

#else
        for (int i = 0; i < ORDER; i++)
        {
            v[i] = x;
            v[i][i]+=delta;
            integrateConst(Osc,v[i],t,t+Osc.T,Osc.dt);
        }
        integrateConst(Osc,x,t,t+Osc.T,Osc.dt,bufferInfo);

#endif // MULTI_THREAD


        const auto endTime = std::chrono::high_resolution_clock::now();

        calTime.push_back({t,duration_cast<duration<double, std::milli>>(endTime-startTime).count()});



        std::array<std::array<double,ORDER>,ORDER> jacobiMatrix;
        for(int i = 0; i < Osc.order; i++)
        {
            for(int j = 0; j < Osc.order; j++)
            {
                jacobiMatrix[i][j] = (v[j][i] - x[i])/delta;
            }
        }

        for(auto& z:spectrum_z)
        {
            state_type dzdt;
            for(size_t i = 0; i < z.size(); i++)
            {
                for(size_t j = 0; j < z.size(); j++)
                {
                    dzdt[i]+=jacobiMatrix[i][j]*z[j];
                }
            }

            z=dzdt;
        }

        for(size_t i = 0; i < spectrum_z.size(); i++)
        {
            for(size_t j = 0; j < i; j++)
            {
                spectrum_z[i] -= spectrum_z[j] * (spectrum_z[i] & spectrum_z[j]);
            }
            lambdaSum[i]+=log(liczWartosc(spectrum_z[i]));
            normVector(spectrum_z[i]);
        }

        if (step % 100 == 0 ) std::cout<<static_cast<int>(t/Osc.T)<<std::endl;

#if TIME_SD
        state_type tempLambda;
        for(size_t i = 0; i < lambdaSum.size(); i++)
        {
            tempLambda[i]=(lambdaSum[i]/(t-Osc.t0));
        }
        lambdaBuffer.push_back(tempLambda);

        if(lambdaBuffer.size() >= 100)
        {
            sd = countSD(lambdaBuffer);
            lambdaBuffer.clear();
            for(auto i:sd) std::cout<<i<<"\t";
            std::cout<<std::endl;

        }
#endif // TIME_SD



        FileL << t;
        for (auto i:lambdaSum) FileL<<";"<<i/(t-Osc.t0);
        FileL << std::endl;
    }

    for (auto i:lambdaSum) std::cout<<i/(t-Osc.t0)<<"\n";


    std::cout<<"Rozpoczynam zapis do pliku\n";

    for (unsigned int i = 0; i < times.size(); i++)
    {
        FileP << times[i];
        for(int j = 0; j < Osc.order; j++)
        {
            FileP << ';' << x_vec[i][j];
        }
        FileP << '\n';
    }

    double sumT=0;
    for(auto i:calTime)
    {
        FileT << i[0] << ';';
        FileT << i[1] << ';';
        FileT << (sumT+=i[1]) << '\n';
    }

    FileP.close();
    FileT.close();
    FileL.close();
}



void biffurAndLambda()
{
   /* std::fstream FileL, FileS;
    FileL.open( "C:\\chwilowy\\prog\\biffurLaypunov.txt", std::ios::out );
    FileS.open( "C:\\chwilowy\\prog\\biffurState.txt", std::ios::out );
    if( FileL.good() == false || FileS.good() == false)
    {
        std::cout<<"blad"<<std::endl;
        //std::cerr << "Error: " << strerror(errno);
        exit(250);
    }

    state_type x(O2ORDER+2);
    inicialValues(x);
    harmOscilator Oscilator;
    std::vector<state_type> x_vec;
    std::vector<double> times;
    std::vector<double> biffVec;
    pushBackBiffurAndNorm bufferInfo( x_vec , times , biffVec, Oscilator);

    int counter = 0;
    for (double biff = BIFF_START; biff < BIFF_STOP; biff+=dBIFF)
    {
        x[O2ORDER] = 0;
        x[O2ORDER+1] = 0;
        x[2]=0.3;
        x[3]=0.1;
        x[4]=0.2;
        x[5]=0.15;
        std::cout<<counter++;
        Oscilator.changeBiffurParametr(biff);

        Oscilator.printInfo();
        double t;
        bool poincareFlag = 0;


        integrateConst(Oscilator , x , S_t0 , S_NouberOfPeriodSkiped*Oscilator.T , Oscilator.dt, Norm(Oscilator));
        double sd1=1, sd2=1;
        state_type lambda1Buffer, lambda2Buffer;
        for(t = S_NouberOfPeriodSkiped*Oscilator.T; sd1 > S_SD_VAL || sd2 > S_SD_VAL || x[O2ORDER]/(t - S_t0)<x[O2ORDER+1]/(t - S_t0); t+=S_NouberOfPeriodCount*Oscilator.T)
        {
            if(poincareFlag)
            {
                integrateConst(Oscilator , x , t , t+S_NouberOfPeriodCount*Oscilator.T , Oscilator.dt, bufferInfo);
            }
            else integrateConst(Oscilator , x , t , t+S_NouberOfPeriodCount*Oscilator.T , Oscilator.dt, Norm(Oscilator));
            lambda1Buffer.push_back(x[O2ORDER]/(t - S_t0));
            lambda2Buffer.push_back(x[O2ORDER+1]/(t - S_t0));

            if(lambda1Buffer.size()>100)
            {
                sd1 = countSD(lambda1Buffer);
                sd2 = countSD(lambda2Buffer);
                if(!poincareFlag && sd1<S_SD_VAL*100)
                {
                    poincareFlag = 1;
                    std::cout<<"\npoicare start write\n";
                }

                std::cout<<"lambda1 = "<<x[O2ORDER]/(t - S_t0)<<"\t\t"<<"sd1 = "<<sd1<<std::endl;
                std::cout<<"lambda2 = "<<x[O2ORDER+1]/(t - S_t0)<<"\t\t"<<"sd2 = "<<sd2<<std::endl;
                std::cout<<std::endl;
                lambda1Buffer.clear();
                lambda2Buffer.clear();
                if(x[O2ORDER]/(t - S_t0)-x[O2ORDER+1]/(t - S_t0)<0.0005 && (sd1 < S_SD_VAL && sd2 < S_SD_VAL)) x[O2ORDER+1]=x[O2ORDER];
            }
        }
        std::cout<<"lambda1 = "<<x[O2ORDER]/(t - S_t0)<<std::endl;
        std::cout<<"lambda2 = "<<x[O2ORDER+1]/(t - S_t0)<<std::endl;


        std::cout<<"Rozpoczynam zapis do pliku\n";

        FileL << biff<<';';
        FileL << x[O2ORDER] / (t - S_t0)<<';';
        FileL << x[O2ORDER+1] / (t - S_t0)<<'\n';
        FileL.flush();



        for (unsigned int i = 0; i < biffVec.size(); i++)
        {
            FileS << biffVec[i]<<';';
            FileS << times[i];
            for(int j = 0; j < O2ORDER; j++)
            {
                FileS << ';' << x_vec[i][j];
            }
            FileS << '\n';
        }
        FileS.flush();
        biffVec.clear();
        times.clear();
        x_vec.clear();
    }

    FileL.close();
    FileS.close();*/
}
