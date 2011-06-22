/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
#include <algorithm>
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
    REAL offset;
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

    REAL kn = par->ktsqr * (par->N->N(par->ktsqr, par->y) - par->offset );

    if (std::abs(ktsqr - par->ktsqr) < 1e-15)
    {
        //cerr << "ktsqr \\approx par->ktsqr and we can't handle this! y=" << par->y
        //   << " par->ktsqr=" << par->ktsqr << " ktsqr: " << ktsqr << " " << LINEINFO << endl;
            /*// Computed analytically
            result += GSL_SIGN(ktsqr - par->ktsqr)*(par->N->N(par->ktsqr, par->y) - exp(par->ktsqr/4.0));
            cout << "result " << result << " instead of " <<
            (ktsqr*par->N->N(ktsqr, par->y) - par->ktsqr*par->N->N(par->ktsqr, par->y))
            / abs(ktsqr - par->ktsqr) << endl;
        ///FIXME: HANDLE NEXT TO LOWEST STEP
        */
    }
    else   
        result += (ktsqr*par->N->N(ktsqr, par->y) - kn)
            / std::abs(ktsqr - par->ktsqr);
    
    result += kn / sqrt( 4.0*SQR(ktsqr) + SQR(par->ktsqr));

    result /= ktsqr;

    if (par->N->RunningCoupling() == MAXK)
        result *= Alphabar_s(std::max(ktsqr, par->ktsqr));
    else if (par->N->RunningCoupling() == MINK)
        result *= Alphabar_s(std::min(ktsqr, par->ktsqr));
    
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
    if (abs(ktsqr - par->ktsqr) < 1e-16)
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

/*
 * Calculate \partial_Y N
 * offset is used in 2nd order Runge Kutta where the derivative
 * \partial_Y N(Y) = F(Y, N(Y)) is evaluated at point
 * F(Y-\delta Y, N(Y) - \delta Y F(Y, N(Y)) )
 * It's default value is 0
 * NB: It can't be used with the kinematical constraint as in that case the
 * equation is not local in Y
 *
 */
REAL BruteForceSolver::RapidityDerivative(REAL ktsqr, REAL y, REAL offset)
{    
    inthelper_bkmom inthelp;
    inthelp.N=this;
    inthelp.y=y;
    inthelp.ktsqr = ktsqr;
    inthelp.offset=offset;
    gsl_function int_helper;

    if (kinematic_constraint==false)
        int_helper.function=&inthelperf_bkmom_noconstraint;
    else
        int_helper.function=&inthelperf_bkmom_constraint;
    int_helper.params=&inthelp;

    REAL result, abserr; 
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(KTSQRINTITERATIONS);

    int status;
    // If ktsqr != ktsqr', we can tell GSL that there is a difficult point
    // at ktsqr'=ktsqr
    if (ktsqr>ktsqrvals[0] and ktsqr<ktsqrvals[ktsqrvals.size()-2] and true==false)
    {
        REAL range[3]; range[0]=ktsqrvals[0];
        range[1]=ktsqr; range[2]=ktsqrvals[ktsqrvals.size()-1];
        status = gsl_integration_qagp(&int_helper, range, 3, 0, KTSQRINTACCURACY,
            KTSQRINTITERATIONS, workspace, &result, &abserr);
    } else
    {
        status=gsl_integration_qag(&int_helper, ktsqrvals[0],
            ktsqrvals[ktsqrvals.size()-2], 0, KTSQRINTACCURACY, 
            KTSQRINTITERATIONS, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    }
    gsl_integration_workspace_free(workspace);
    if (status ) cerr << "Error " << status << " at " << LINEINFO << ":"
        << " ktsqr=" << ktsqr <<", y=" << y << " result=" << result << ", abserror=" <<
        abserr << " relerror: " << abserr/result << endl;

    // Nonlinear term
    result -= SQR(N(ktsqr, y));

    if (RunningCoupling()==CONSTANT)
        result*=0.2;
    if (RunningCoupling()==PARENT_DIPOLE or RunningCoupling()==MAXK
        or RunningCoupling() == MINK)
        result *= Alphabar_s(ktsqr);
    
    return result;
}
/*
 * GSL ODEIV evolution
 * Returns RapidityDerivative for every ktsqr in vector
 * TODO: Doesn't probably work with kinematical constraint, does it?
 *
 */
struct EvolutionHelper
{
    BruteForceSolver* N;
};
int Evolve(REAL y, const REAL amplitude[], REAL result[], void *params)
{
    cout << "Evolving with y=" << y << endl;
    EvolutionHelper* helper = (EvolutionHelper*)params;
   // Actually we don't need amplitude[] as it only contains
   // elements N[ktsqrind][yind], but GSL uses it to determine whether
   // the volution is accurate 
   #pragma omp parallel for
   for( int i=0 ; i<helper->N->KtsqrPoints() ; i++) {
      result[i] = helper->N->RapidityDerivative(helper->N->Ktsqrval(i), y);
   }
   return GSL_SUCCESS;
}


// Solve BK, lowest order
void BruteForceSolver::Solve(REAL maxy)
{
    if (averages>0)
    {
        cerr << "Averagement is disabled..." << endl;
    }
    // Find maxyind corresponding to maxy
    int maxyind=YPoints();
    for (unsigned int i=1; i<=YPoints(); i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }
    int largedifference=0;

    // **** used in GSL solver *****
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;

    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, KtsqrPoints());
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0, 0.05);    //abserr relerr
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (KtsqrPoints());
    EvolutionHelper help; help.N=this;
    gsl_odeiv_system sys = {Evolve, NULL, KtsqrPoints(), &help};
    REAL Y=0;
    REAL Yi=0;
    REAL h = 0.1;   // Original step size
    REAL *amplitude=new REAL[KtsqrPoints()];
    for (int ktsqrind=0; ktsqrind < KtsqrPoints(); ktsqrind++)
    {
        amplitude[ktsqrind] = n[ktsqrind][0];
    }
    
    // ******************************
    for (int yind=1; yind<=maxyind; yind++)
    {
        cerr << "Solving for y=" << yvals[yind] << endl;
        // Solve N(y+DELTA_Y, kt) for every kt

        // Use RungeKutta = GSL ODEIV system, doesn't work (yet?) with kin.
        // constraint
        if (rungekutta)
        {

            Yi = yvals[yind];
            // gsl_odeiv_evolve_apply increases Y according to the step size
            // why gsl_odeiv_evolve_apply doesn't set Y=Y_i at the end? 
            while (Y<Yi)
            {
                int status = gsl_odeiv_evolve_apply (e, c, s, &sys, 
                                &Y, Yi, &h, amplitude);
                if (status != GSL_SUCCESS) {
                    cerr << "Error in gsl_odeiv_evolve_apply at " << LINEINFO
                        << ": " << gsl_strerror(status) << " (" << status << ")"
                        << " y=" << Y << ", h=" << h << endl;
             }
                cout << "Evolved up to " << Y << "/" << Yi << ", h=" << h << endl;
            }
            cout << "Solved yind " << yind << " to Y=" << Y << " with step size "
                << h << endl;
            
      
            if (std::abs(Y-Yi)>0.05)
            {
                cerr << "Y-Y_i is " << Y-Yi << ", Y=" << Y << ", Yi=" << Yi << endl;
            }
            for (int ktsqrind=0; ktsqrind < KtsqrPoints(); ktsqrind++)
            {
                n[ktsqrind][yind] = amplitude[ktsqrind];
            }

            

        }
        else
        {


            #pragma omp parallel for
            for (int ktsqrind=0; ktsqrind<KtsqrPoints(); ktsqrind++)
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

                if (rungekutta and yind>1)
                {
                    /*REAL der2 = RapidityDerivative(tmpkt, yvals[yind-1] - dy,
                            - dy*tmpder);
                    newn = n[ktsqrind][yind-1]
                        + 3.0/2.0*dy*tmpder - 1.0/2.0*dy*der2;
                    */
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

        /*
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
                for (int ktsqrind=0; ktsqrind<KtsqrPoints()-1; ktsqrind++)
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
        } */
    
    }

    delete[] amplitude;

}

BruteForceSolver::BruteForceSolver()
{
    rungekutta = false;
}
BruteForceSolver::~BruteForceSolver()
{

}

void BruteForceSolver::SetRungeKutta(bool rk)
{
    rungekutta = rk;
}
