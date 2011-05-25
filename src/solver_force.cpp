/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
#include "solver_force.hpp"
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <cmath>
using std::abs;

// Amplitude class methods which solve BK directly (=brute force)


/*
 * \partial_Y from BK in momentum space
 */

struct inthelper_bkmom
{
    REAL ktsqr;
    Amplitude* N;
    REAL y;
};


/*
 * BK equation in momentum space integrated over \theta (b-indep. situation)
 * Ref e.g. hep-ph/0110325 eq (8)
 */
REAL inthelperf_bkmom_noconstraint(REAL ktsqr, void* p)
{
    inthelper_bkmom* par = (inthelper_bkmom*) p;

    REAL result=0;

    /* When ktsqr -> par->ktsqr first doesn't diverge, but numerics would fail
     *  as there is term 1/(ktsqr - par->ktsqr)
     * Analytically one can expand N(k')=N(k+\eps) and find that in that limit
     * the integrand is
     *  sgn(k'^2 - k^2)/k^2 [N(k^2) + k'^2 N'(k^2) ] + N(k^2)/( Sqrt[5] k^2 )
     */

    if (std::abs(ktsqr - par->ktsqr) < 1e-14)
    {
        cerr << "ktsqr \\approx par->ktsqr and we can't handle this! y=" << par->y
           << " par->ktsqr=" << par->ktsqr << " ktsqr: " << ktsqr << endl;
            /*// Computed analytically
            result += GSL_SIGN(ktsqr - par->ktsqr)*(par->N->N(par->ktsqr, par->y) - exp(par->ktsqr/4.0));
            cout << "result " << result << " instead of " <<
            (ktsqr*par->N->N(ktsqr, par->y) - par->ktsqr*par->N->N(par->ktsqr, par->y))
            / abs(ktsqr - par->ktsqr) << endl;
        ///FIXME: HANDLE NEXT TO LOWEST STEP
        */
    }
    else   
        result += (ktsqr*par->N->N(ktsqr, par->y) - par->ktsqr*par->N->N(par->ktsqr, par->y))
            / std::abs(ktsqr - par->ktsqr);
    
    result += par->ktsqr * par->N->N(par->ktsqr, par->y) / sqrt( 4.0*SQR(ktsqr) + SQR(par->ktsqr));

    result /= ktsqr;
    
    //cout << ktsqr << "  " << result << endl;
    return result;
}

/*
 * BK equation in momentum space integrated over \theta (b-indep. situation)
 * With kinematical constraint z < k^2/k'^2, or ln 1/z > ln k'^2/k^2
 * => Y-y > ln k'^2/k^2 when integrating over y up to Y
 * Ref e.g. hep-ph/0110325 eq (44)
 * Differentiating that w.r.t. y gives
 * \partial_Y N(k,Y) = \bar \alpha_s \int d^2 k'/k'^2 {
 *     *[ \theta(k^2-k'^2)k'^2 N(k',Y) + \theta(k'^2-k^2) k'^2 N(k', Y-ln(k'^2/k^2)) - k^2 N(k,y)]
 *     / Abs[k^2-k'^2]
 *     + k^2 N(k,Y) / Sqrt[4k'^4+k^4] } - \bar \alpha_s N(k,Y)^2
 * This computes the integrand
 */
REAL inthelperf_bkmom_constraint(REAL ktsqr, void* p)
{
    inthelper_bkmom* par = (inthelper_bkmom*) p;

    REAL result=0;
    if (abs(ktsqr - par->ktsqr) < 1e-14)
    {
        cerr << "ktsqr \\approx par->ktsqr and we can't handle this! y=" << par->y
            << " par->ktsqr=" << par->ktsqr << " at " << LINEINFO << endl;
            ///TODO
    }
    else
    {
        if (par->ktsqr > ktsqr)
            result += ktsqr*par->N->N(ktsqr, par->y);
        else
            result += ktsqr*par->N->N(ktsqr, par->y - log(ktsqr/par->ktsqr));
        result -= par->ktsqr*par->N->N(par->ktsqr, par->y);
        result /= abs(ktsqr - par->ktsqr);
    }
    result += par->ktsqr*par->N->N(par->ktsqr, par->y)/sqrt(4.0*SQR(ktsqr)+SQR(par->ktsqr));

    result/=ktsqr;

    return result;
}

REAL BruteForceSolver::RapidityDerivative(REAL ktsqr, REAL y)
{
    REAL alphabar = 0.2;    //TODO
    inthelper_bkmom inthelp;
    inthelp.N=this;
    inthelp.y=y;
    inthelp.ktsqr = ktsqr;
    gsl_function int_helper;

    if (kinematic_constraint==false)
        int_helper.function=&inthelperf_bkmom_noconstraint;
    else
        int_helper.function=&inthelperf_bkmom_constraint;
    int_helper.params=&inthelp;

    REAL result, abserr; 
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(KTSQRINTITERATIONS); 
    int status=gsl_integration_qag(&int_helper, ktsqrvals[0], ktsqrvals[ktsqrvals.size()-2], 0, KTSQRINTACCURACY, 
        KTSQRINTITERATIONS, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    if (status or result>1e6 or result<-1e6) cerr << "Error " << status << " at " << LINEINFO << ":"
        << " ktsqr=" << ktsqr <<", y=" << y << " result=" << result << ", abserror=" <<
        abserr << " relerror: " << abserr/result << endl;

    // Nonlinear term
    result -= SQR(N(ktsqr, y));
    
    return alphabar*result;
}


// Solve BK, lowest order
void BruteForceSolver::Solve(REAL maxy)
{
    // Find maxyind corresponding to maxy
    int maxyind=YPoints();
    for (unsigned int i=1; i<=YPoints(); i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }
    int largedifference=0;
    for (int yind=1; yind<=maxyind; yind++)
    {
        cerr << "Solving for y=" << yvals[yind] << endl;
        // Solve N(y+DELTA_Y, kt) for every kt

#pragma omp parallel for
        for (unsigned int ktsqrind=0; ktsqrind<KtsqrPoints()-1; ktsqrind++)
        {
            REAL tmpkt = ktsqrvals[ktsqrind];
            REAL dy = yvals[yind]-yvals[yind-1];
            REAL tmpder = RapidityDerivative(tmpkt, yvals[yind-1]);
            REAL newn=n[ktsqrind][yind-1] + dy*tmpder;

            // Adams method: apprximate derivative as a 2nd order polynomial
            // Y_{i+1} = y_i + hf_i + 1/2 h (f_i - f_{i-1} )
            // We can use this only for yind>1
            if (yind>1 and adams_method==true)
            {
                REAL old_der = derivatives[ktsqrind][yind-2];
                newn = n[ktsqrind][yind-1] + dy*tmpder
                    + 1.0/2.0*dy*( tmpder - old_der);
            }
            
            
            AddDataPoint(ktsqrind, yind, newn, tmpder );
            if( abs(newn - n[ktsqrind][yind-1])/n[ktsqrind][yind-1] > 0.1)
            {
                largedifference++;
            }
            //cout << "N(ktsqr=" << ktsqrvals[ktsqrind] <<", y=" << yvals[yind] << ") = " << newn
            //<< ",  at lower rapidity it was " << n[ktsqrind][yind-1] << endl;
            
        }
    }
    cout << endl << "#" << largedifference << " out of " << maxyind * (KtsqrPoints()-1)
            << " too large differences" << endl;
    // Again
    for (int avg=0; avg<averages; avg++)
    {
        
        largedifference=0;
        for (int yind=1; yind<maxyind-1; yind++)
        {
            //cout << "Solving for y=" << yvals[yind] << endl;
            // Solve N(y+DELTA_Y, kt) for every kt
            ///TODO: Different iterations are not independent, but the difference
            /// caused by different order of execution should be higher order?
#pragma omp parallel for
            for (unsigned int ktsqrind=0; ktsqrind<KtsqrPoints()-1; ktsqrind++)
            {
                REAL tmpkt = ktsqrvals[ktsqrind];
                REAL tmpder = RapidityDerivative(tmpkt, yvals[yind]);
                REAL dy = yvals[yind]-yvals[yind-1];
                REAL newn = n[ktsqrind][yind-1] + dy*0.5*(tmpder+derivatives[ktsqrind][yind-1]);

                if( abs(newn - n[ktsqrind][yind-1])/n[ktsqrind][yind-1] > 0.02)
                {
                    largedifference++;
                }

                AddDataPoint(ktsqrind, yind, newn, 0.5*(tmpder + derivatives[ktsqrind][yind-1]));

                
            }
        }
        cout << endl << "#" << largedifference << " out of " << maxyind * (KtsqrPoints()-1)
            << " too large differences" << endl;
    }

}

BruteForceSolver::BruteForceSolver()
{

}
BruteForceSolver::~BruteForceSolver()
{

}
