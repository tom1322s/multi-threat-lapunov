#ifndef HARMOSC_H
#define HARMOSC_H

#include <array>
#define ORDER 6

#define MULTI_THREAD 1

//time - 0 sdCount -1
#define TIME_SD 1

#define OSC_KIND 0

// ³apac x3 kiedy x1 jest 0
// ³apac okres z ile razy przejdzie przez 0 i pomijac okresy przejsciowe

typedef std::array<double,ORDER> state_type;

#if OSC_KIND == 0
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
        const double sdValRef = 1e-7;
        double sdVal = 1e-7;
        const int nouberOfPeriodSkipedRef = 300.0;
        int nouberOfPeriodSkiped = 300.0;
        const int nouberOfPeriod = 4000;
        const int nouberOfPeriodCount = 10;
        //nouberOfPeriodPoincare 200;
        const double t0 = 0.0;


        const double BIFF_START = 0.0;
        const double BIFF_STOP = 4.0;
        const double dBIFF = ((BIFF_STOP - BIFF_START) / 100.0);

        //const double dtVal=37890.0/10.0; // !!!!!!!!

    protected:

    private:

        int biff_type;
        double k,d,a;
};
#elif OSC_KIND == 1
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
        const double sdValRef = 1e-8;
        double sdVal = 1e-8;
        const int nouberOfPeriodSkipedRef = 300.0;
        int nouberOfPeriodSkiped = 300.0;
        const int nouberOfPeriod = 4000;
        const int nouberOfPeriodCount = 10;
        //nouberOfPeriodPoincare 200;
        const double t0 = 0.0;


        const double BIFF_START = 1.1;
        const double BIFF_STOP = 1.4;
        const double dBIFF = ((BIFF_STOP - BIFF_START) / 100.0);

        //const double dtVal=37890.0/10.0; // !!!!!!!!

    protected:

    private:

        int biff_type;
        double h,q,eta;
};

#else
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
        const double sdValRef = 1e-7;
        double sdVal = 1e-7;
        const int nouberOfPeriodSkipedRef = 300.0;
        int nouberOfPeriodSkiped = 300.0;
        const int nouberOfPeriod = 4000;
        const int nouberOfPeriodCount = 10;
        //nouberOfPeriodPoincare 200;
        const double t0 = 0.0;


        const double BIFF_START = 0.94;
        const double BIFF_STOP = 0.97;
        const double dBIFF = ((BIFF_STOP - BIFF_START) / 5.0);

        //const double dtVal=37890.0/10.0; // !!!!!!!!

    protected:

    private:

        int biff_type;
        double ksiY, q, eta, ni, delta, ksi;
};

#endif // OSC_KIND

#endif // HARMOSC_H
