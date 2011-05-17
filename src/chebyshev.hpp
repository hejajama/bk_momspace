/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _CHEBYSHEV_H
#define _CHEBYSHEV_H

/*
 * Chebyshev polynomial based basis
 * This class represents a vector for which the components are coefficients of
 * the Chebyshev polynomials
 * In that space we define an inner product
 * a \cdot b = \int_{-1}^1 (1-x)^{-1/2} a(x)b(x)
 */

#include "config.hpp"
#include <gsl/gsl_chebyshev.h>
#include <vector>

class ChebyshevBasis
{
    public:
        ChebyshevBasis(unsigned int d);
        REAL Component(unsigned int n);
        REAL InnerProduct(ChebyshevBasis &vec);
        REAL Evaluate(REAL x);
        void SetComponent(unsigned int c, REAL val);


    private:
        unsigned int degree;
        gsl_cheb_series *cheb;

};

#endif
