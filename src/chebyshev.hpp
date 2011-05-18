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
#include <iostream>

class ChebyshevVector
{
    public:
        ChebyshevVector(unsigned int d);
        ChebyshevVector(std::vector<REAL> vec);
        void SetLimits(REAL a, REAL b);
        REAL Component(unsigned int n);
        REAL DotProduct(ChebyshevVector &vec);
        // Dot product with an arbitrary function
        REAL DotProduct( REAL(*f)(REAL x, void* p), void* p );
        REAL Evaluate(REAL x);
        void SetComponent(unsigned int c, REAL val);
        void Normalize();

        ChebyshevVector operator+(ChebyshevVector &v);
        ChebyshevVector operator-(ChebyshevVector &v);
        ChebyshevVector& operator*(REAL x);
        ChebyshevVector& operator=(ChebyshevVector v);

        unsigned int Degree();
    private:
        unsigned int degree;
        gsl_cheb_series *cheb;

};

std::ostream& operator<<(std::ostream& os, ChebyshevVector& v);

#endif
