/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "spectrum.hpp"
#include "amplitude.hpp"
#include <tools/tools.hpp>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

const int MAXITER_KTINT = 2000;
const int MAXITER_THETAINT = 200;
const REAL INTACCURACY_KT = 0.01;

/* dN_ch / dydp_T, integrated over b and angular dep. of p_T
 * Requires normalization info from other source
 *
 * Ref: arXiv:0707:2545: (integrated over angular dep. of p_T and imp. par.
 * dN_ch / dy dp_t
 * = C 1/p_t^2 \int_0^p_t dk_T k_T d\theta \alpha_s(Q)
 *   * \psi(x1, 0.5|k_t+p_T| ) * \psi(x2, 0.5*|k_T - p_T| )
 * Q = 0.5 max{k_t \pm p_t}
 * x_{1,2} = p_T/sqrt(s) exp(\pm y)
 *
 * Then we make change of variables to u = ln kt (NB! kt, not ktsqr)
 */
struct Inthelper_kt
{
    REAL lnkt;
    REAL theta;
    Amplitude* N;
    Spectrum* Spect;
    REAL y; REAL sqrts; REAL pt;
};
REAL Inthelperf_kt_k(REAL kt, void* p);
REAL Inthelperf_kt_theta(REAL theta, void* p);
REAL Spectrum::dNch_dydpsqr(REAL sqrts, REAL y, REAL psqr)
{
    REAL x1 = std::sqrt(psqr)/sqrts*std::exp(y);
    REAL x2 = std::sqrt(psqr)/sqrts*std::exp(-y);

    REAL y1 = N->Y(x1);
    REAL y2 = N->Y(x2);
    //cout << "p_T=" << std::sqrt(psqr) << " y=" << y <<
    //      " eval_y1=" << N->Y(std::sqrt(psqr)/sqrts*exp(y)) << " x1=" << x1 << 
    //      " eval_y2=" << N->Y(std::sqrt(psqr)/sqrts*exp(-y)) << " x2=" << x2 << endl;
    
    REAL p = std::sqrt(psqr);
    Inthelper_kt helper;
    helper.N = N;
    helper.y = y; helper.sqrts=sqrts; helper.pt=p;

    gsl_function int_helper;
    int_helper.function=&Inthelperf_kt_k;
    int_helper.params=&helper;

    gsl_integration_workspace *workspace 
            = gsl_integration_workspace_alloc(MAXITER_KTINT);
    REAL minlnk = std::log(std::sqrt(N->Ktsqrval(0)));
    REAL maxlnk = std::log(p);
    REAL result, abserr;
    int status = gsl_integration_qag(&int_helper, minlnk,
        maxlnk, 0, INTACCURACY_KT, MAXITER_KTINT, GSL_INTEG_GAUSS41, workspace,
        &result, &abserr);
    if (status)
    {
        cerr << "kt integral failed at " << LINEINFO
            <<", pt = " << helper.pt << ", result " << result << " relerr "
            << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
    return result;
    
}
// Integral helper functions
REAL Inthelperf_kt_k(REAL lnkt, void *p)
{
    Inthelper_kt *par = (Inthelper_kt*) p;
    par->lnkt = lnkt;

    gsl_function int_helper;
    int_helper.function = &Inthelperf_kt_theta;
    int_helper.params=p;
    gsl_integration_workspace *workspace 
            = gsl_integration_workspace_alloc(MAXITER_THETAINT);
    REAL mintheta = 0;
    REAL maxtheta = 2.0*M_PI;
    REAL result, abserr;
    int status = gsl_integration_qag(&int_helper, mintheta, maxtheta,
        0, INTACCURACY_KT, MAXITER_THETAINT, GSL_INTEG_GAUSS41, workspace,
        &result, &abserr);
    gsl_integration_workspace_free(workspace);
    if (status)
    {
        cerr << "theta integral failed at " << LINEINFO << ": lnkt=" << lnkt
            <<", pt = " << par->pt << ", result " << result << " relerr "
            << std::abs(abserr/result) << endl;
    }
    return result;
}
REAL Inthelperf_kt_theta(REAL theta, void* p)
{
    Inthelper_kt* par = (Inthelper_kt*) p;
    REAL pt = par->pt;
    REAL kt = std::exp(par->lnkt);
    REAL result=0;
    // |k_T + p_T|/2
    REAL plus = std::sqrt( SQR(pt) + SQR(kt) + 2.0*kt*pt*std::cos(theta) ) / 2.0;
    // |k_T - p_T|/2
    REAL minus = std::sqrt( SQR(pt) + SQR(kt) - 2.0*kt*pt*std::cos(theta) ) / 2.0;

    REAL Q = std::max(plus, minus);

    REAL x1 = pt/par->sqrts*std::exp(par->y);
    REAL x2 = pt/par->sqrts*std::exp(-par->y);

    REAL y1 = par->N->Y(x1);
    REAL y2 = par->N->Y(x2);
    
    return 1.0/SQR(pt)*kt*kt*Alpha_s(SQR(Q))
        *par->N->N(plus, y1)*par->N->N(minus, y2)*std::pow(1.0-x1,4)*std::pow(1.0-x2,4);    
}

/*
 * Charged hadron multiplicity as a function of y
 * = Spectrum::dNch_dydpsqr integrated over p_T
 * Angular integral gives just 2\pi, but we have dropped the normalization
 * constants everywhere
 */

REAL Inthelperf_pt(REAL pt, void* p);

REAL Spectrum::dNch_dy(REAL sqrts, REAL y)
{
    Inthelper_kt help;
    help.N=N; help.y=y; help.sqrts=sqrts;
    help.Spect=this;

    gsl_function int_helper;
    int_helper.function=&Inthelperf_pt;
    int_helper.params=&help;

    gsl_integration_workspace *workspace 
            = gsl_integration_workspace_alloc(MAXITER_KTINT);
    REAL minlnp = std::log(std::sqrt(N->Ktsqrval(0)));
    REAL maxlnp = std::max(std::log(std::sqrt(N->Ktsqrval(N->KtsqrPoints()-1))),
                            std::log(sqrts/2.0) );
    REAL result, abserr;
    int status = gsl_integration_qag(&int_helper, minlnp,
        maxlnp, 0, INTACCURACY_KT, MAXITER_KTINT, GSL_INTEG_GAUSS41, workspace,
        &result, &abserr);
    if (status)
    {
        cerr << "pt integral failed at " << LINEINFO
            <<", y = " << y << ", result " << result << " relerr "
            << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
    return result;

}

REAL Inthelperf_pt(REAL lnpt, void* p)
{
    Inthelper_kt* par = (Inthelper_kt*) p;
    // pt -> lnpt => dpt pt -> dlnpt exp(2lnpt)
    return std::exp(2.0*lnpt)*par->Spect->dNch_dydpsqr(par->sqrts, par->y,
        SQR(std::exp(lnpt)) );
}

Spectrum::Spectrum(Amplitude* amp)
{
    N = amp;
}

