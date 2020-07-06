#ifndef HARMOSC_H
#define HARMOSC_H

#include <vector>
typedef std::vector< double > state_type;

class HarmOsc
{
    public:
        HarmOsc();
        void printInfo();
        bool checkPoincareRequirement(const double t, const double tOld);
        double calculateTimeToSavePoint(const double t, const double tOld);
        void changeBiffurParametr(double biffVal);
        double getBiffurParametr();
        void updatePeriod();

        void operator() (const state_type &x , state_type &dxdt , const double  t  );

        double T, dt;

        const int order = 6;
        const double sdVal = 1e-9;
        const int nouberOfPeriodSkiped = 400;//300.0
        const int nouberOfPeriodCount = 10;
        //nouberOfPeriodPoincare 200;
        const double t0 = 0.0;


        const double BIFF_START = 0.0;
        const double BIFF_STOP = 4.0;
        const double dBIFF = ((BIFF_STOP - BIFF_START) / 100.0);

        const double dtVal=37890.0/10.0; // !!!!!!!!

    protected:

    private:

        int biff_type;
        double k,d,a;
};

#endif // HARMOSC_H
