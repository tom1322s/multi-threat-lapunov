#ifndef OBSERVER_H_INCLUDED
#define OBSERVER_H_INCLUDED

#include <vector>
#include <array>
//typedef std::vector< double > state_type;

//double liczWartosc(state_type &x, int numberPer);
//void normVector(state_type &x, int numberPer);
//void normPerturbationAndUpdatePerpendicularComponent(state_type &x);
//void copyPerturbation(state_type &vec1, state_type &vec2, const state_type &scr);
//state_type operator * (state_type& scr, double scalar);
//state_type& operator -= (state_type &vec, const state_type &scr);
//double operator & (const state_type &vec1, const state_type &vec2);

/*
struct Obs
{
    std::vector< double >& m_times;
    std::vector< state_type >& m_states;
    std::vector< double > & m_lambda;
    std::vector< double > & m_lambda2;
    double tTemp;
    int coutCounter;
    harmOscilator &Osc;



    Obs( std::vector< double > &times, std::vector< state_type > &states, std::vector< double > &lambda , std::vector< double > &lambda2, harmOscilator &Osc )
    : m_times( times ), m_states( states ) , m_lambda( lambda ) , m_lambda2( lambda2 ) , tTemp(0) , coutCounter(1), Osc(Osc){ }

    void operator()( state_type &x , double t )
    {
        if((t - tTemp) > Osc.T/100)
        {
            m_times.push_back( t );
            m_states.push_back( x );
            m_lambda.push_back(x[O2ORDER]/(t - S_t0));
            m_lambda2.push_back(x[O2ORDER+1]/(t - S_t0));
            tTemp = t;
        }
        //Osc.updateJacobian(x,t);
        normPerturbationAndUpdatePerpendicularComponent(x);

        if(t > coutCounter * Osc.T)
        {
            std::cout<<coutCounter<<std::endl;
            coutCounter++;
        }
    }
};*/
template<typename OSC>
struct pushBackTimeAndState
{
    std::vector< double >& times;
    std::vector< std::vector< double > >& states;
    double tTemp;
    OSC &Osc;

    pushBackTimeAndState( std::vector< double > &times, std::vector< std::vector< double > > &states, OSC &Osc )
    : times( times ), states( states ) , tTemp(0) , Osc(Osc){ }

    template <typename T>
    void operator()( T &x , double t )
    {
        if((t - tTemp) > Osc.T/100)
        {
            times.push_back( t );
            std::vector<double> temp;
            temp.reserve(x.size());
            for(auto i:x) temp.push_back(i);
            states.push_back( temp );
            tTemp = t;
        }
    }
};


/*
struct Norm
{
    harmOscilator &Osc;
    double tTemp;
    int counter;

    Norm(harmOscilator &Osc)
   :  Osc( Osc ),tTemp(0), counter(1) { }

    void operator()( state_type &x , double t )
    {
        normPerturbationAndUpdatePerpendicularComponent(x);

        if((t - tTemp) > Osc.T && t<S_NouberOfPeriodSkiped*Osc.T)
        {
            std::cout<<counter++<<std::endl;
            //std::cout<<t<<std::endl;
            tTemp = t;
        }
    }
};*/

/*struct pushBackBiffurAndNorm
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;
    std::vector< double >& m_biff;
    harmOscilator& m_Oscilator;
    double tOld;
    state_type xOld;
    double tTemp;
    int coutCounter;


    pushBackBiffurAndNorm ( std::vector< state_type > &states , std::vector< double > &times ,std::vector< double >& biff, harmOscilator& Oscilator )
    : m_states( states ) , m_times( times ) , m_biff( biff ) , m_Oscilator( Oscilator ) , tTemp(0) , coutCounter(1){ }

    void operator()( state_type &x , double t )
    {
        if( m_Oscilator.checkPoincareRequirement(t, tOld))
        {
            double tSave = m_Oscilator.calculateTimeToSavePoint(t, tOld);
            state_type xSave;
            for (int i = 0; i < O2ORDER; i++)
            {
                xSave.push_back(x[i] - ((t - tSave) * (x[i] - xOld[i])) / (t - tOld));
            }
            m_states.push_back( xSave );
            m_times.push_back( tSave );
            m_biff.push_back( m_Oscilator.getBiffurParametr());
        }
        normPerturbationAndUpdatePerpendicularComponent(x);
        tOld = t;
        xOld = x;

       // if((t - tTemp) > m_Oscilator.T)
       // {
       //     std::cout<<coutCounter++<<std::endl;
       //     //std::cout<<t<<std::endl;
       //     tTemp = t;
       // }

    }
};

*/

#endif // OBSERVER_H_INCLUDED
