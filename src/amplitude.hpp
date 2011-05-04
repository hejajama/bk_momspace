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

const REAL DELTA_Y=0.1;
const int MAXY=50;
const REAL KTSQR_MULTIPLIER = 1.02;  // ktsqr_{i+1} = ktsqr_i*KTSQR_MULTIPLIER
const REAL MINKTSQR = 1e-7;
const REAL MAXKTSQR = 1e9;

const int POINTS_Y= (int)(MAXY/DELTA_Y);
const int POINTS_KTSQR = (int)(std::log(MAXKTSQR/MINKTSQR) / std::log(KTSQR_MULTIPLIER) );

const REAL KTSQRINTACCURACY = 0.01;
const int KTSQRINTITERATIONS = 7000;

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
        REAL N(REAL ktsqr, REAL y);
        REAL RapidityDerivative(REAL ktsqr, REAL y);

        void AddDataPoint(int ktsqrindex, int yindex, REAL value, REAL der);

        void Interpolate();

        void Solve(REAL maxy);

        REAL Ktsqrval(int i);
        REAL Yval(int i);
        void SetInitialCondition(INITIAL_CONDITION i);
        void SetKinematicConstraint(bool kc);
        string InitialConditionStr();
        
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

       // doublexyz *interpolated_amplitude;

        
        
        
};

#endif
