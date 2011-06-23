/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _SOLVER_FORCE_H
#define _SOLVER_FORCE_H

#include "config.hpp"
#include "amplitude.hpp"
#include <gsl/gsl_integration.h>
#include <vector>
#include <cmath>

class BruteForceSolver : public Amplitude
{
    public:
        void Solve(REAL maxy);
        REAL RapidityDerivative(REAL ktsqr, REAL y);
        
        void SetRungeKutta(bool rk);
        void SetAdamsMethod(bool ad);

        BruteForceSolver();
        ~BruteForceSolver();
    private:
        bool rungekutta;

        void InitializeAdamsMethod();

        bool adams_method;   // Whether the Adams method should be used when solvin DE


};

#endif
