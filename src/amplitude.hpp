/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _AMPLITUDE_H
#define _AMPLITUDE_H

#include "config.hpp"
#include <vector>
#include <cmath>
//#include <bci.h>

const REAL DEFAULT_DELTA_Y=0.1;
const unsigned int DEFAULT_MAXY=50;
const REAL DEFAULT_KTSQR_MULTIPLIER = 1.015;  // ktsqr_{i+1} = ktsqr_i*KTSQR_MULTIPLIER
const REAL DEFAULT_MINKTSQR = 1e-8;
const REAL DEFAULT_MAXKTSQR = 1e10;

//const unsigned int POINTS_Y= (int)(MAXY/DELTA_Y);
//const unsigned int POINTS_KTSQR = (int)(std::log(MAXKTSQR/MINKTSQR) / std::log(KTSQR_MULTIPLIER) );

const REAL KTSQRINTACCURACY = 0.005;
const int KTSQRINTITERATIONS = 9000;
const int INTERPOLATION_POINTS = 10;

enum INITIAL_CONDITION
{
    FTIPSAT,    // \int d^2r 1/(2\pi r^2) exp(ik.r) 2(1-exp(-r^2))
    INVPOWER    // 1/(k^2 + 1), as in BK in full mom. space, hep-ph/0504080
};
    

class Amplitude
{
    public:
        Amplitude();
        void Initialize();
		void Clear();
        REAL N(REAL ktsqr, REAL y);
        REAL RapidityDerivative(REAL ktsqr, REAL y);
        REAL LogLogDerivative(REAL ktsqr, REAL y);

        void AddDataPoint(int ktsqrindex, int yindex, REAL value, REAL der);

        int ReadData(string file);

        void Interpolate();

        void Solve(REAL maxy);

        REAL Ktsqrval(unsigned int i);
        REAL Yval(unsigned int i);
        void SetInitialCondition(INITIAL_CONDITION i);
        void SetKinematicConstraint(bool kc);
        string InitialConditionStr();
        void SetNumberOfAveragements(int avg);

        void SetMinKtsqr(REAL mkt);
        void SetMaxKtsqr(REAL mkt);
        void SetKtsqrMultiplier(REAL m);
        void SetMaxY(REAL y);
        void SetDeltaY(REAL dy);

        unsigned int YPoints();
        unsigned int KtsqrPoints();
        REAL KtsqrMultiplier();
        REAL DeltaY();
        
    private:
        REAL InitialCondition(REAL ktsqr);  // N() at y=0

        // Values of N as a function of ktsqr and y
        std::vector<REAL> ktsqrvals;
        std::vector<REAL> yvals;
        // n[ktsqr][y]
        std::vector< std::vector<REAL> > n;

        // derivatives
        std::vector< std::vector<REAL> > derivatives;

        INITIAL_CONDITION ic;
        bool kinematic_constraint;
        int averages;

        REAL minktsqr;
        REAL maxktsqr;
        REAL ktsqr_multiplier;
        REAL maxy;
        REAL delta_y;

        bool datafile;	// True if data is read from external file

        bool adams_method;   // Whether the Adams method should be used when solvin DE
        
        
};

#endif
