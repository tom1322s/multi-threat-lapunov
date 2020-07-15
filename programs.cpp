#include "programs.h"
#include "functions.h"
#include <thread>
#include <future>
#include <fstream>
//#include <math.h>
#include <chrono>
#include <set>

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
    //std::vector<double> times;
    //std::vector<state_type> x_vec;
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
        std::set<double> x_P;
        pushBackPoicare bufferInfo(x_P);
        state_type lambdaSum = {};
        std::cout<<counter++;
        Osc.changeBiffurParametr(biff);

        Osc.printInfo();
        double t=0.0;

        //
        if (liczWartosc(x)<1e-3)
        {
            inicialValues(x);
        }

        /*for(auto& i:spectrum_z)
        {
            inicialValues(i);
        }*/
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

#if TIME_SD
            if (sd[0]<Osc.sdVal*10)
#else
            if (step > 0.9*Osc.nouberOfPeriod)
#endif // TIME_SD
            {
                integrateConst(Osc,x,t,t+Osc.T,Osc.dt,bufferInfo);
            }
            else
            {
                integrateConst(Osc,x,t,t+Osc.T,Osc.dt);
            }

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
#if TIME_SD
            if (sd[0]<Osc.sdVal*10)
#else
            if (step > 0.9*Osc.nouberOfPeriod)
#endif // TIME_SD
            {
                integrateConst(Osc,x,t,t+Osc.T,Osc.dt,bufferInfo);
            }
            else
            {
                integrateConst(Osc,x,t,t+Osc.T,Osc.dt);
            }

#endif // MULTI_THREAD



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

        for (auto point:x_P)
        {
            FileP << biff;
            FileP << ';' << point;
            FileP << '\n';
        }

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

void timePhraseLaypunov2()
{
    std::fstream FileL, FileT, FileP;
    FileP.open( "C:\\chwilowy\\prog\\2timePhrase.txt", std::ios::out );
    FileL.open( "C:\\chwilowy\\prog\\2timeLaponov.txt", std::ios::out );
    FileT.open( "C:\\chwilowy\\prog\\2diff.txt", std::ios::out );
    if( FileP.good() == false ||  FileT.good() == false || FileL.good() == false)
    {
        std::cout<<"blad"<<std::endl;
        //std::cerr << "Error: " << strerror(errno);
        exit(250);
    }

    FileT<<"dx0;kwx0";
    for(int i =1;i<ORDER;i++) FileT<<";dx"<<i<<";kwx"<<i;
    for(int i =0;i<ORDER;i++)
    {
        FileT<<";dL"<<i<<";"<<"kwL"<<i;
    }
    FileT<<std::endl;

    state_type xSumK;
    for(auto& i:xSumK) i=0.0;
    state_type lamSumK;
    for(auto& i:lamSumK) i=0.0;

    HarmOsc Osc;
    std::vector<double> times1;
    std::vector<std::vector<double>> x_vec1;
    pushBackTimeAndState<HarmOsc> bufferInfo1( times1, x_vec1,Osc);
    std::vector<double> times2;
    std::vector<std::vector<double>> x_vec2;
    pushBackTimeAndState<HarmOsc> bufferInfo2( times2, x_vec2,Osc);

    state_type x1;
    inicialValues(x1);
    state_type x2;
    x2=x1;

    std::array<state_type,ORDER> spectrum_z1;
    for(auto& i:spectrum_z1)
    {
        inicialValues(i);
    }
    std::array<state_type,ORDER> spectrum_z2;
    spectrum_z2=spectrum_z1;

    Osc.printInfo();
    double delta = 1e-6;
    state_type lambdaSum1;
    state_type lambdaSum2;
    std::vector<std::vector<double>> diffPar;

    double t1=0.0;
    double t2=0.0;

#if TIME_SD
    state_type sd1;
    for(auto& i:sd1) i = 1.0;
    std::vector<state_type> lambdaBuffer1;
    state_type sd2;
    for(auto& i:sd2) i = 1.0;
    std::vector<state_type> lambdaBuffer2;
#endif // TIME_SD

#if TIME_SD
    for (long int step = 0; step < Osc.nouberOfPeriodSkiped || isSdGood(sd1, Osc.sdVal); step++)
#else
    for (long int step = 0; step < Osc.nouberOfPeriod; step++)
#endif // TIME_SD
    {
        t1 = Osc.t0 + step * Osc.T;


        std::array<state_type,ORDER> v1;
        std::array<state_type,ORDER> v2;


        std::array<std::future<state_type>,ORDER> thread1;
        std::array<std::future<state_type>,ORDER> thread2;
        for (int i = 0; i < Osc.order; i++)
        {
            v1[i]=x1;
            v1[i][i]+=delta;
            thread1[i] = std::async(std::launch::async, &integrateConstR<HarmOsc,state_type>, Osc,v1[i],t1,t1+Osc.T,Osc.dt);
            v2[i]=x2;
            v2[i][i]+=delta;
            thread2[i] = std::async(std::launch::async, &integrateConstRStare<HarmOsc,state_type>, Osc,v2[i],t2,t2+Osc.T,Osc.dt);
        }
        integrateConst(Osc,x1,t1,t1+Osc.T,Osc.dt,bufferInfo1);
        integrateConstStare(Osc,x2,t2,t2+Osc.T,Osc.dt,bufferInfo2);
        for (int i = 0; i < Osc.order; i++)
        {
            v1[i] = thread1[i].get();
        }
        for (int i = 0; i < Osc.order; i++)
        {
            v2[i] = thread2[i].get();
        }




        std::array<std::array<double,ORDER>,ORDER> jacobiMatrix1;
        for(int i = 0; i < Osc.order; i++)
        {
            for(int j = 0; j < Osc.order; j++)
            {
                jacobiMatrix1[i][j] = (v1[j][i] - x1[i])/delta;
            }
        }

        for(auto& z:spectrum_z1)
        {
            state_type dzdt;
            for(size_t i = 0; i < z.size(); i++)
            {
                for(size_t j = 0; j < z.size(); j++)
                {
                    dzdt[i]+=jacobiMatrix1[i][j]*z[j];
                }
            }

            z=dzdt;
        }

        for(size_t i = 0; i < spectrum_z1.size(); i++)
        {
            for(size_t j = 0; j < i; j++)
            {
                spectrum_z1[i] -= spectrum_z1[j] * (spectrum_z1[i] & spectrum_z1[j]);
            }
            lambdaSum1[i]+=log(liczWartosc(spectrum_z1[i]));
            normVector(spectrum_z1[i]);
        }

        std::array<std::array<double,ORDER>,ORDER> jacobiMatrix2;
        for(int i = 0; i < Osc.order; i++)
        {
            for(int j = 0; j < Osc.order; j++)
            {
                jacobiMatrix2[i][j] = (v2[j][i] - x2[i])/delta;
            }
        }

        for(auto& z:spectrum_z2)
        {
            state_type dzdt;
            for(size_t i = 0; i < z.size(); i++)
            {
                for(size_t j = 0; j < z.size(); j++)
                {
                    dzdt[i]+=jacobiMatrix2[i][j]*z[j];
                }
            }

            z=dzdt;
        }

        for(size_t i = 0; i < spectrum_z2.size(); i++)
        {
            for(size_t j = 0; j < i; j++)
            {
                spectrum_z2[i] -= spectrum_z2[j] * (spectrum_z2[i] & spectrum_z2[j]);
            }
            lambdaSum2[i]+=log(liczWartosc(spectrum_z2[i]));
            normVector(spectrum_z2[i]);
        }

        if (step % 100 == 0 ) std::cout<<static_cast<int>(t1/Osc.T)<<std::endl;

#if TIME_SD
        state_type tempLambda;
        for(size_t i = 0; i < lambdaSum1.size(); i++)
        {
            tempLambda[i]=(lambdaSum1[i]/(t1-Osc.t0));
        }
        lambdaBuffer1.push_back(tempLambda);

        if(lambdaBuffer1.size() >= 100)
        {
            sd1 = countSD(lambdaBuffer1);
            lambdaBuffer1.clear();
            //for(auto i:sd) std::cout<<i<<"\t";
            //std::cout<<std::endl;

        }

        for(size_t i = 0; i < lambdaSum2.size(); i++)
        {
            tempLambda[i]=(lambdaSum2[i]/(t2-Osc.t0));
        }
        lambdaBuffer2.push_back(tempLambda);

        if(lambdaBuffer2.size() >= 100)
        {
            sd2 = countSD(lambdaBuffer2);
            lambdaBuffer2.clear();
            //for(auto i:sd) std::cout<<i<<"\t";
            //std::cout<<std::endl;

        }
#endif // TIME_SD



        FileL << t1;
        for (auto i:lambdaSum1) FileL<<";"<<i/(t1-Osc.t0);
        FileL << t2;
        for (auto i:lambdaSum2) FileL<<";"<<i/(t2-Osc.t0);
        FileL << std::endl;



        for(int i = 0; i<ORDER; i++)
        {
            xSumK[i]+=pow(x2[i]-x1[i],2);
            lamSumK[i]+=pow((lambdaSum2[i]/(t2-Osc.t0))-(lambdaSum1[i]/(t1-Osc.t0)),2);
        }

        FileT<<x2[0]-x1[0]<<";"<<xSumK[0];
        for(int i = 1; i<ORDER; i++) FileT<<";"<<x2[i]-x1[i]<<";"<<xSumK[i];
        for(int i = 0; i<ORDER; i++)
        {
            FileT<<";"<<(lambdaSum2[i]/(t2-Osc.t0))-(lambdaSum1[i]/(t1-Osc.t0))<<";"<<xSumK[i];
        }
        FileT << std::endl;

        t2 += Osc.T;
    }

    //for (auto i:lambdaSum) std::cout<<i/(t-Osc.t0)<<"\n";


    std::cout<<"Rozpoczynam zapis do pliku\n";

    for (unsigned int i = 0; i < times1.size(); i++)
    {
        FileP << times1[i];
        for(int j = 0; j < Osc.order; j++)
        {
            FileP << ';' << x_vec1[i][j];
        }
        FileP << times2[i];
        for(int j = 0; j < Osc.order; j++)
        {
            FileP << ';' << x_vec2[i][j];
        }
        FileP << '\n';
    }



    FileP.close();
    FileT.close();
    FileL.close();
}

