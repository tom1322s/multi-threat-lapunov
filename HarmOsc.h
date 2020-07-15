#ifndef HARMOSC_H
#define HARMOSC_H

#include <array>
#define ORDER 6

#define MULTI_THREAD 1

//time - 0 sdCount -1
#define TIME_SD 1

// ³apac x3 kiedy x1 jest 0
// ³apac okres z ile razy przejdzie przez 0 i pomijac okresy przejsciowe

typedef std::array<double,ORDER> state_type;

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

        const int order = ORDER;
        const double sdVal = 1e-7;
        const int nouberOfPeriodSkiped = 300.0;
        const int nouberOfPeriod = 1000;
        const int nouberOfPeriodCount = 10;
        //nouberOfPeriodPoincare 200;
        const double t0 = 0.0;


        const double BIFF_START = 0.0;
        const double BIFF_STOP = 4.0;
        const double dBIFF = ((BIFF_STOP - BIFF_START) / 1000.0);

        //const double dtVal=37890.0/10.0; // !!!!!!!!

    protected:

    private:

        int biff_type;
        double k,d,a;
};

#endif // HARMOSC_H
