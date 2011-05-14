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

const int MAXITER_VINT=1000;
const int MAXITER_UINT=MAXITER_VINT;
const REAL UVINTACCURACY=0.001;



struct Integrand_helper
{
    ChebyshevAmplitudeSolver* N;
    REAL v;
    REAL y;
    REAL n;
    REAL m;
};

REAL Integrand_helperu(REAL v, void* p);
REAL Integrand_helperv(REAL v, void *p)
{
    Integrand_helper* h = (Integrand_helper*) p;
    h->v = v;
    
    gsl_function int_helper;
    int_helper.function=&Integrand_helperu;
    int_helper.params=h;
    REAL result, abserr;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MAXITER_UINT);
    int status=gsl_integration_qag(&int_helper, -1, 1, 0, UVINTACCURACY, 
        MAXITER_UINT, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (status) std::cerr << "Error " << status << " at " << LINEINFO
        << ": Result " << result << ", abserror: " << abserr 
        << " (v=" << v <<")" << endl;
    return result;
    
    
}

REAL Integrand_helperu(REAL u, void* p)
{
    Integrand_helper* h = (Integrand_helper*) p;
    REAL result=0;
    REAL n=h->n;
    REAL m=h->m;
    REAL v=h->v;
    if (std::abs(u-v)<eps) return 0.0;  //TODO: CHECK
    result = h->N->Chebyshev(n, v) - std::exp(h->N->Delta(u,v))*h->N->Chebyshev(n, u);
    result /= std::abs(exp(h->N->Delta(u,v))-1);
    result += h->N->Chebyshev(n, u)/std::sqrt(1.0+4.0*std::exp(-2.0*h->N->Delta(u,v)));
    result /= std::sqrt(1.0-SQR(u));
    
    return result;
    
}

void ChebyshevAmplitudeSolver::Solve(REAL maxy)
{
    Prepare();

    /* Compute matrix F_{nm}
     * = \int_0^1 dv \int_-1^1 du (1-u^2)^(-1/2) T_m(u) {
     * (T_n(v) - exp(-delta)T_n(u) ) / Abs(exp(delta)-1)
     * + T_n(u)/Sqrt(1+4*exp(-2*delta)) }
     * here delta = (m1+m2)(u-v)
     */

    Integrand_helper h;
    h.n=1;
    h.m=2;
    h.N=this;
    h.y=0;

    gsl_function int_helper;
    int_helper.function=&Integrand_helperv;
    int_helper.params=&h;
    REAL result, abserr;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MAXITER_UINT);
    int status=gsl_integration_qag(&int_helper, 0, 1, 0, UVINTACCURACY, 
        MAXITER_UINT, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (status) std::cerr << "Error " << status << " at " << LINEINFO
        << ": Result " << result << ", abserror: " << abserr << endl;

    cout << result << " pm " << abserr << endl;



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
        REAL* tmpc = new REAL[CHEBYSHEV_DEGREE];
        coef.push_back(tmpc);
    }

    gsl_cheb_series *cs = gsl_cheb_alloc (CHEBYSHEV_DEGREE);
    ICHelper help;
    help.N=this;
    gsl_function f;
    f.function = ICHelperf;
    f.params=&help;

    gsl_cheb_init(cs, &f, 0.0, 1.0);

    for (unsigned int i=0; i<CHEBYSHEV_DEGREE; i++)
        coef[0][i]=cs->c[i];

    gsl_cheb_free(cs);

}


ChebyshevAmplitudeSolver::~ChebyshevAmplitudeSolver()
{
    for (unsigned int yind=0; yind<YPoints(); yind++)
    {
        delete[] coef[yind];
    }
    gsl_cheb_free(cheb);

}

ChebyshevAmplitudeSolver::ChebyshevAmplitudeSolver()
{
    cheb = gsl_cheb_alloc (CHEBYSHEV_DEGREE);
    cheb->a=-1.0; cheb->b=1.0;
    for (unsigned int i=1; i<=CHEBYSHEV_DEGREE; i++)
        cheb->c[i]=0;
    cheb->c[0]=2.0;
    
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

/*
 * Evaluate T_n(x)
 */
REAL ChebyshevAmplitudeSolver::Chebyshev(unsigned int n, REAL x)
{
    if (n > CHEBYSHEV_DEGREE) 
    {
        cerr << "Asked Chebyshev polynomial at degree n=" << n 
            << ", maximum  is " << CHEBYSHEV_DEGREE << endl;
        return -1;
    }
    if (oldn==n)
        return gsl_cheb_eval(cheb, x);
    
    cheb->c[oldn]=0.0;
    if (n==0)
    {
        cheb->c[0]=2.0;
    }
    else 
    {
        cheb->c[n]=1.0;
    }
    
    REAL result = gsl_cheb_eval_n(cheb, n, x);
    oldn=n;
    return result;
    
    
}

REAL ChebyshevAmplitudeSolver::Delta(REAL u, REAL v)
{
    return (m1+m2)*(u-v);    
}
