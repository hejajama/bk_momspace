/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _SOLVER_FORCE2_H
#define _SOLVER_FORCE2_H

#include "config.hpp"
#include "amplitude.hpp"
#include <gsl/gsl_integration.h>
#include <vector>
#include <cmath>

// Same as BruteForceSolver but uses compactified real axis

class BruteForceSolver2 : public Amplitude
{
    public:
        void Solve(REAL maxy);
        REAL UAmplitude(REAL u, REAL y);
        REAL RapidityDerivative(REAL u, REAL y);

        REAL Ktsqr(REAL u);
        REAL U(REAL ktsqr);
        
        

        BruteForceSolver2();
        ~BruteForceSolver2();
    private:
        void Prepare();
    
        std::vector< std::vector<REAL> > amplitude; // n[u][y]

        int upoints;
        REAL minu, maxu;
        std::vector<REAL> uvals;


};

#endif
