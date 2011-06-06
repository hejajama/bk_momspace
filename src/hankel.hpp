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
#include <gsl/gsl_dht.h>

class Hankel
{
    public:
        Hankel(Amplitude *amp, REAL y_, int npoints);
        void PrintRAmplitude();  // Print amplitude in position space
        ~Hankel();
    private:
        int points;
        Amplitude* N;
        REAL *sample;
        REAL *transformed;
        gsl_dht* dht;
        REAL y;

    
};

#endif
