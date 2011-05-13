/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
#include "chebyshev_amplitude.hpp"
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_chebyshev.h>

#include <cmath>
using std::abs;


void ChebyshevAmplitudeSolver::Solve(REAL maxy)
{
    Prepare();

    /* Compute matrix F_{nm}
     * = \int_0^1 dv \int_-1^1 du (1-u^2)^(-1/2) T_m(u) {
     * (T_n(v) - exp(-delta)T_n(u) ) / Abs(exp(delta)-1)
     * + T_n(u)/Sqrt(1+4*exp(-2*delta)) }
     * here delta = (m1+m2)(u-v)
     */


    



};

/*
 * Prepare everything
 * Allocate memory and compute the coefficients for the Chebyshev approx
 * of the initial condition
 */

// Helper structures we use to compute Chebyshev expansion of the
// initial condition
struct ICHelper
{
    ChebyshevAmplitudeSolver* N;
};

REAL ICHelperf(REAL x, void* p)
{
    ICHelper* par = (ICHelper*) p;
    REAL ktsqr = std::exp( ( par->N->M1() + par->N->M2() )*x - par->N->M1() );
    return par->N->InitialCondition(ktsqr);   
}

void ChebyshevAmplitudeSolver::Prepare()
{
    m1 = -std::log(MinKtsqr());
    m2 = std::log(MaxKtsqr());
    for (unsigned int yind=0; yind < YPoints(); yind++)
    {
        REAL* tmpc = new REAL[CHEBYSHEV_ORDER];
        coef.push_back(tmpc);
    }

    gsl_cheb_series *cs = gsl_cheb_alloc (CHEBYSHEV_ORDER);
    ICHelper help;
    help.N=this;
    gsl_function f;
    f.function = ICHelperf;
    f.params=&help;

    gsl_cheb_init(cs, &f, 0.0, 1.0);

    for (unsigned int i=0; i<CHEBYSHEV_ORDER; i++)
        coef[0][i]=cs->c[i];

    gsl_cheb_free(cs);

}


ChebyshevAmplitudeSolver::~ChebyshevAmplitudeSolver()
{
    for (unsigned int yind=0; yind<YPoints(); yind++)
    {
        delete[] coef[yind];
    }

}

REAL ChebyshevAmplitudeSolver::M1()
{
    return m1;
}

REAL ChebyshevAmplitudeSolver::M2()
{
    return m2;
}

REAL ChebyshevAmplitudeSolver::Ktsqr(REAL u)
{
    // Comput ktsqr from u
    return std::exp( (m1+m2)*u - m1 );
}
