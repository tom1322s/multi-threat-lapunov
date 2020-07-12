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
    //std::vector<std::vector<double>> lambda;
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
    std::fstream FileL, FileT, FileP;
    FileP.open( "C:\\chwilowy\\prog\\bifurPoicare.txt", std::ios::out );
    FileL.open( "C:\\chwilowy\\prog\\bifurLaponov.txt", std::ios::out );
    FileT.open( "C:\\chwilowy\\prog\\bifurtime.txt", std::ios::out );
    if( FileP.good() == false ||  FileT.good() == false || FileL.good() == false)
    {
        std::cout<<"blad"<<std::endl;
        //std::cerr << "Error: " << strerror(errno);
        exit(250);
    }

    HarmOsc Osc;
    std::vector<double> times;
    std::vector<state_type> x_vec;
    //std::vector<std::vector<double>> lambda;

    state_type x;
    inicialValues(x);

    std::array<state_type,ORDER> spectrum_z;
    for(auto& i:spectrum_z)
    {
        inicialValues(i);
    }

    double delta = 1e-6;
    std::vector<std::vector<double>> calTime;


    int counter = 0;
    for (double biff = Osc.BIFF_START; biff < Osc.BIFF_STOP; biff+=Osc.dBIFF)
    {
        state_type lambdaSum = {};
        std::cout<<counter++;
        Osc.changeBiffurParametr(biff);

        Osc.printInfo();
        double t=0.0;

        //
        inicialValues(x);
        for(auto& i:spectrum_z)
        {
            inicialValues(i);
        }
        //

        #if TIME_SD
        state_type sd;
        for(auto& i:sd) i = 1.0;
        std::vector<state_type> lambdaBuffer;
#endif // TIME_SD

        const auto startTime = std::chrono::high_resolution_clock::now();

#if TIME_SD
        for (long int step = 0; step < Osc.nouberOfPeriodSkiped || isSdGood(sd, Osc.sdVal); step++)
#else
        for (long int step = 0; step < Osc.nouberOfPeriod; step++)
#endif // TIME_SD
        {
            t = Osc.t0 + step * Osc.T;
            std::array<state_type,ORDER> v;

#if MULTI_THREAD

            std::array<std::future<state_type>,ORDER> thread;
            for (int i = 0; i < Osc.order; i++)
            {
                v[i]=x;
                v[i][i]+=delta;
                thread[i] = std::async(std::launch::async, &integrateConstR<HarmOsc,state_type>, Osc,v[i],t,t+Osc.T,Osc.dt);
            }
            integrateConst(Osc,x,t,t+Osc.T,Osc.dt);
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
            integrateConst(Osc,x,t,t+Osc.T,Osc.dt);

#endif // MULTI_THREAD

#if TIME_SD
            //if (step > 0.8*Osc.nouberOfPeriodSkiped)
            if (sd[0]<Osc.sdVal*10)
#else
            if (step > 0.9*Osc.nouberOfPeriod)
#endif // TIME_SD
            {
                times.push_back(t);
                x_vec.push_back(x);
            }



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
                //for(auto i:sd) std::cout<<i<<"\t";
                //std::cout<<std::endl;

            }
#endif // TIME_SD


        }
        const auto endTime = std::chrono::high_resolution_clock::now();

        calTime.push_back({biff,duration_cast<duration<double, std::milli>>(endTime-startTime).count()});

        FileL << biff;
        for (auto i:lambdaSum) FileL<<";"<<i/(t-Osc.t0);
        FileL << std::endl;

        for (unsigned int i = 0; i < times.size(); i++)
        {
            FileP << biff;
            FileP << ';' << times[i];
            for(int j = 0; j < Osc.order; j++)
            {
                FileP << ';' << x_vec[i][j];
            }
            FileP << '\n';
        }
        times.clear();
        x_vec.clear();

    }

    std::cout<<"Rozpoczynam zapis do pliku\n";



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
