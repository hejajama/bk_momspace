/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

// Calculate discrete Hankel transformation using GSL
// So get N(r) from N(k) via
// N(r) = r^2 \int_0^\infty k*J_0(k*r)*N(k)

#ifndef _HANKEL_HPP
#define _HANKEL_HPP

#include "config.hpp"
#include "amplitude.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sum.h>
#include <cmath>

class Hankel
{
    public:
        Hankel(Amplitude *amp);
        REAL Amplitude_r(REAL r, REAL y);
        ~Hankel();
    private:
        int points;
        Amplitude* N;
        REAL *sample;
        REAL *transformed;
        static const double epsilon1=1.0e-12;
        static const double epsilon=1.0e-12;
        REAL y;

    
};

#endif
