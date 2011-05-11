/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef SOLVER_CHEBYSHEV_H
#define SOLVER_CHEBYSHEV_H

#include "amplitude.hpp"
#include "config.hpp"
#include <vector>
using std::vector;
/*
 * Solves BK equation by writing the solution in the basis of
 * Chebyshev polynomials
 */

const unsigned int CHEBYSHEV_DEGREE = 40;


class ChebyshevSolver : public Amplitude
{

    public:
        void Solve(REAL maxy);

    private:
        // \int_0^1 dv f(v) using Chebyshev polynomial approximation
        // v is the parameter integrated over, u is fixed
        REAL Integrate(int uind, int yind);
        REAL Integrand(int uind, int vind, int yind); 
        // \int_0^1 dx T_n(x)
        REAL ChebyshevIntegral(int n);
        void Prepare();
        vector<REAL> uvals;       // u grid, u \in [0,1]

        // L = ln k^2
        REAL minl, maxl;




};


#endif
