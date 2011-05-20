/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _CHEBYSHEV_AMPLITUDE_H
#define _CHEBYSHEV_AMPLITUDE_H

#include "config.hpp"
#include "amplitude.hpp"
#include "chebyshev.hpp"
#include <vector>
#include <gsl/gsl_chebyshev.h>
#include <cmath>

enum BASIS_BOUNDARY_CONDITION
{
    CHEBYSHEV,
    CHEBYSHEV_ZERO,        // Chebyshevs which are 0 at x=1
};


class ChebyshevAmplitudeSolver : public Amplitude
{
    public:
        ChebyshevAmplitudeSolver();
        ~ChebyshevAmplitudeSolver();
        void Solve(REAL maxy);
        void Prepare();
        REAL Chebyshev(unsigned int n, REAL x);
        REAL Delta(REAL u, REAL v);
        REAL Basis(unsigned int n, REAL x);     // Evaluate ith basis function
        ChebyshevVector& BasisVector(unsigned int n); // Poitner to nth basis vec

        REAL N(REAL ktsqr, REAL y);

        REAL M1();
        REAL M2();
        REAL Ktsqr(REAL u);
        REAL U(REAL ktsqr);

        void SetChebyshevDegree(unsigned int d);
        void SetBoundaryCondition(BASIS_BOUNDARY_CONDITION bc);
        


    private:
        // Chebyshev polynomial coefficients, vector index is yind
        // Each pointer points to an array with CHEBYSHEV_ORDER numbers
        // This for so it is efficient to use with GSL
        std::vector< std::vector<REAL> > coef;

        // Basis vectors
        std::vector< ChebyshevVector > basis;
        void ComputeBasisVectors();
        
        REAL m1, m2;    // ln k^2=L, L \in [-M1, M2]
        
        unsigned int oldn;
        gsl_cheb_series *cheb;    // Table used in Chebyshev evaluation function
        gsl_integration_qaws_table* qaws_table;
        unsigned int chebyshev_degree;
        static const unsigned int DEFAULT_CHEBYSHEV_DEGREE=40;
        BASIS_BOUNDARY_CONDITION boundary_condition;

};

// Helper structures we use to compute Chebyshev expansion of the
// initial condition
struct ICHelper
{
    ChebyshevAmplitudeSolver* N;
};


#endif
