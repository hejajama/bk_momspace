/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _CHEBYSHEV_AMPLITUDE_H
#define _CHEBYSHEV_AMPLITUDE_H

#include "config.hpp"
#include "amplitude.hpp"
#include <vector>
#include <gsl/gsl_chebyshev.h>
#include <cmath>



class ChebyshevAmplitudeSolver : public Amplitude
{
    public:
        ChebyshevAmplitudeSolver();
        ~ChebyshevAmplitudeSolver();
        void Solve(REAL maxy);
        void Prepare();
        REAL Chebyshev(unsigned int n, REAL x);
        REAL Delta(REAL u, REAL v);

        REAL M1();
        REAL M2();
        


    private:
        // Chebyshev polynomial coefficients, vector index is yind
        // Each pointer points to an array with CHEBYSHEV_ORDER numbers
        // This for so it is efficient to use with GSL
        std::vector< REAL* > coef;
        REAL m1, m2;    // ln k^2=L, L \in [-M1, M2]
        REAL Ktsqr(REAL u);
        unsigned int oldn;
        gsl_cheb_series *cheb;    // Table used in Chebyshev evaluation function
        static const unsigned int CHEBYSHEV_DEGREE=40;

};


#endif
