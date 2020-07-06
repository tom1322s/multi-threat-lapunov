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
    std::fstream File;
    File.open( "C:\\chwilowy\\prog\\timePhrase.txt", std::ios::out );
    if( File.good() == false )
    {
        std::cout<<"blad"<<std::endl;
        //std::cerr << "Error: " << strerror(errno);
        exit(250);
    }

    HarmOsc Osc;
    std::vector<double> times;
    std::vector<state_type> x_vec;
    pushBackTimeAndState<HarmOsc> bufferInfo( times, x_vec,Osc);

    state_type x(Osc.order);
    inicialValues(x);

    std::vector<state_type> spectrum_z;
    for (int i = 0; i < Osc.order; i++)
    {
        state_type z(Osc.order);
        inicialValues(z);
        spectrum_z.push_back(std::move(z));
    }

    Osc.printInfo();
    double delta = 1e-6;
    state_type lambdaSum(Osc.order);

    for (double t = Osc.t0; t < Osc.nouberOfPeriodSkiped*Osc.T; t+=Osc.T)
    {
        std::vector<state_type> v;
        std::future<state_type> thread[Osc.order];
        for (int i = 0; i < Osc.order; i++)
        {
            v.push_back(state_type (x));
            v[i][i]+=delta;
            thread[i] = std::async(std::launch::async, &integrateConst<HarmOsc>, Osc,v[i],t,t+Osc.T,Osc.dt);
        }
        integrateConst(Osc,x,t,t+Osc.T,Osc.dt,bufferInfo);
        for (int i = 0; i < Osc.order; i++)
        {
            v[i] = thread[i].get();
        }

       // for (int i = 0; i < ORDER; i++)
       // {
       //     v.push_back(state_type (x));
       //     v[i][i]+=delta;
       //     integrateConst(Oscilator,v[i],t,t+Oscilator.T,Oscilator.dt,bufferInfo);
       // }
       // integrateConst(Oscilator,x,t,t+Oscilator.T,Oscilator.dt,bufferInfo);



        std::array<std::array<double,2>,2> jacobiMatrix; // !!!!!!!!!!!!!
        for(int i = 0; i < Osc.order; i++)
        {
            for(int j = 0; j < Osc.order; j++)
            {
                jacobiMatrix[i][j] = (v[j][i] - x[i])/delta;
            }
        }

        /*for(int i = 0; i < ORDER; i++)
        {
            for(int j = 0; j < ORDER; j++)
            {
                std::cout<<jacobiMatrix[i][j]<<"\t";
            }
            std::cout<<std::endl;
        }*/

        //for_each(z.begin(),z.end(),[](double i){std::cout<<i<<"\t";});
        //std::cout<<std::endl;


        for(auto& z:spectrum_z)
        {

            state_type dzdt(z.size());
            for(int i = 0; i < z.size(); i++)
            {
                for(int j = 0; j < z.size(); j++)
                {
                    dzdt[i]+=jacobiMatrix[i][j]*z[j];
                }
            }

            z=std::move(dzdt);
        }

        for(int i = 0; i < spectrum_z.size(); i++)
        {
            for(int j = 0; j < i; j++)
            {
                spectrum_z[i] -= spectrum_z[j] * (spectrum_z[i] & spectrum_z[j]);
            }
            lambdaSum[i]+=log(liczWartosc(spectrum_z[i]));
            normVector(spectrum_z[i]);
        }


        std::cout<<static_cast<int>(t/Osc.T)<<std::endl;
        for (auto i:lambdaSum) std::cout<<i/t<<"\t\t";
        std::cout<<std::endl;

    }

    /*double sd1=1, sd2=1;
    state_type lambda1Buffer, lambda2Buffer;
    for(double t = S_NouberOfPeriodSkiped*Oscilator.T; sd1 > S_SD_VAL || sd2 > S_SD_VAL || x[O2ORDER]/(t - S_t0)<=x[O2ORDER+1]/(t - S_t0); t+=S_NouberOfPeriodCount*Oscilator.T)
    {
        integrateConst( Oscilator , x , t , t+S_NouberOfPeriodCount*Oscilator.T , Oscilator.dt, bufferInfo);
        lambda1Buffer.push_back(x[O2ORDER]/(t - S_t0));
        lambda2Buffer.push_back(x[O2ORDER+1]/(t - S_t0));

        if(lambda1Buffer.size()>100)
        {
            sd1 = countSD(lambda1Buffer);
            sd2 = countSD(lambda2Buffer);

            std::cout<<"lambda1 = "<<x[O2ORDER]/(t - S_t0)<<"\t\t"<<"sd1 = "<<sd1<<std::endl;
            std::cout<<"lambda2 = "<<x[O2ORDER+1]/(t - S_t0)<<"\t\t"<<"sd2 = "<<sd2<<std::endl;
            std::cout<<std::endl;
            lambda1Buffer.clear();
            lambda2Buffer.clear();
            //std::cout<<"#############################"<<std::endl;
        }
        if(x[O2ORDER]/(t - S_t0)-x[O2ORDER+1]/(t - S_t0)<0.0005 && (sd1 < S_SD_VAL && sd2 < S_SD_VAL)) x[O2ORDER+1]=x[O2ORDER];
    }*/

    std::cout<<"Rozpoczynam zapis do pliku\n";

    for (unsigned int i = 0; i < times.size(); i++)
    {
        File << times[i];
        for(int j = 0; j < Osc.order; j++)
        {
            File << ';' << x_vec[i][j];
        }
        //File << ';' << lambda[i];
        //File << ';' << lambda2[i];
        File << '\n';
    }

    File.close();
}

void timePhraseLaypunovSDCount()
{
    std::fstream File;
    File.open( "C:\\chwilowy\\prog\\timePhrase.txt", std::ios::out );
    if( File.good() == false )
    {
        std::cout<<"blad"<<std::endl;
        //std::cerr << "Error: " << strerror(errno);
        exit(250);
    }

    HarmOsc Osc;
    std::vector<double> times;
    std::vector<state_type> x_vec;
    pushBackTimeAndState<HarmOsc> bufferInfo( times, x_vec,Osc);

    state_type x(Osc.order);
    inicialValues(x);

    std::vector<state_type> spectrum_z;
    for (int i = 0; i < Osc.order; i++)
    {
        state_type z(Osc.order);
        inicialValues(z);
        spectrum_z.push_back(std::move(z));
    }

    Osc.printInfo();
    double delta = 1e-6;
    state_type lambdaSum(Osc.order);

    std::vector<double> sd(Osc.order,1.0);
    std::vector<std::vector<double>> lambdaBuffer(Osc.order);

    for (double t = Osc.t0; t < Osc.nouberOfPeriodSkiped*Osc.T || isSdGood(sd, Osc.sdVal); t+=Osc.T)
    {
        //const auto startTime = std::chrono::high_resolution_clock::now();

        //t=Osc.t0;
        //std::vector<state_type> v;
        //std::future<state_type> thread[Osc.order];
        //std::vector<std::future<state_type>> thread(Osc.order);
        //for (int i = 0; i < Osc.order; i++)
        //{
        //    v.push_back(state_type (x));
        //    v[i][i]+=delta;
        //    thread[i] = std::async(std::launch::async, &integrateConst<HarmOsc>, Osc,v[i],t,t+Osc.T,Osc.dt);
       // }
        //auto Osc2=Osc;
        //integrateConst(Osc,x,t,t+Osc.T,Osc.dt,bufferInfo);
        //x = integrateConst(Osc,x,t,t+Osc.T,Osc.dt);
        auto mainThread = std::async(std::launch::async, &integrateConst<HarmOsc>, Osc,x,t,t+Osc.T,Osc.dt);
        //for (int i = 0; i < Osc.order; i++)
       // {
       //     v[i] = thread[i].get();
        //}
        x = mainThread.get();

        //const auto endTime = std::chrono::high_resolution_clock::now();

        //std::cout<<(duration_cast<duration<double, std::milli>>(endTime-startTime).count());



      /*  for (int i = 0; i < Osc.order; i++)
        {
            v.push_back(state_type (x));
            v[i][i]+=delta;
            integrateConst(Osc,v[i],t,t+Osc.T,Osc.dt,bufferInfo);
        }
        integrateConst(Osc,x,t,t+Osc.T,Osc.dt,bufferInfo);*/

/*
        std::vector<std::vector<double>> jacobiMatrix(x.size());
        for(auto& i:jacobiMatrix)
        {
            i.resize(x.size());
        }

        for(int i = 0; i < Osc.order; i++)
        {
            for(int j = 0; j < Osc.order; j++)
            {
                jacobiMatrix[i][j] = (v[j][i] - x[i])/delta;
            }
        }

                //for(int i = 0; i < ORDER; i++)
                //{
                //   for(int j = 0; j < ORDER; j++)
                //    {
                //        std::cout<<jacobiMatrix[i][j]<<"\t";
                //    }
                //    std::cout<<std::endl;
                //}

        //for_each(z.begin(),z.end(),[](double i){std::cout<<i<<"\t";});
        //std::cout<<std::endl;


        for(auto& z:spectrum_z)
        {

            state_type dzdt(z.size());
            for(int i = 0; i < z.size(); i++)
            {
                for(int j = 0; j < z.size(); j++)
                {
                    dzdt[i]+=jacobiMatrix[i][j]*z[j];
                }
            }

            z=std::move(dzdt);
        }

        for(int i = 0; i < spectrum_z.size(); i++)
        {
            for(int j = 0; j < i; j++)
            {
                spectrum_z[i] -= spectrum_z[j] * (spectrum_z[i] & spectrum_z[j]);
            }
            lambdaSum[i]+=log(liczWartosc(spectrum_z[i]));
            normVector(spectrum_z[i]);
        }

*/
        //static int counter=0;
       // if (counter==10)
       // {

       // counter = 0;
        //std::cout<<static_cast<int>(t/Osc.T)<<std::endl;
        //}
       // else counter++;

   /*     for (auto i:lambdaSum) std::cout<<i/(t-Osc.t0)<<"\t\t";
        std::cout<<std::endl;

        for(int i = 0; i < lambdaSum.size(); i++)
        {
            lambdaBuffer[i].push_back(lambdaSum[i]/(t-Osc.t0));
        }

        if(lambdaBuffer[0].size() >= 50)
        {
            std::cout<<"#################"<<std::endl<<std::endl;
            for (int i = 0; i < lambdaBuffer.size(); i++)
            {
                sd[i] = countSD(lambdaBuffer[i]);
                std::cout<<sd[i]<<"\t";
                lambdaBuffer[i].clear();
            }
            std::cout<<std::endl<<std::endl<<"#################"<<std::endl;

        }
*/

    }

    std::cout<<"Rozpoczynam zapis do pliku\n";

    for (unsigned int i = 0; i < times.size(); i++)
    {
        File << times[i];
        for(int j = 0; j < Osc.order; j++)
        {
            File << ';' << x_vec[i][j];
        }
        //File << ';' << lambda[i];
        //File << ';' << lambda2[i];
        File << '\n';
    }

    File.close();
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
