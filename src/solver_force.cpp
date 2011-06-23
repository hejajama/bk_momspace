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
#include <gsl/gsl_odeiv2.h> // Requires GSL 1.15
#include <cmath>
using std::abs;

// Amplitude class methods which solve BK directly (=brute force)


/*
 * \partial_Y from BK in momentum space
 */

struct inthelper_bkmom
{
    REAL lnktsqr;
    Amplitude* N;
    REAL y;
    REAL offset;
};


/*
 * BK equation in momentum space integrated over \theta (b-indep. situation)
 * Ref e.g. hep-ph/0110325 eq (8)
 */
REAL inthelperf_bkmom_noconstraint(REAL lnktsqr, void* p)
{
    inthelper_bkmom* par = (inthelper_bkmom*) p;

    REAL result=0;
    REAL ktsqr = std::exp(lnktsqr);
    REAL parktsqr = std::exp(par->lnktsqr);

    /* When lnktsqr -> par->lnktsqr first doesn't diverge, but numerics would fail
     *  as there is term 1/(ktsqr - par->ktsqr)
     * Analytically one can expand N(k')=N(k+\eps) and find that in that limit
     * the integrand is
     *  sgn(k'^2 - k^2)/k^2 [N(k^2) + k'^2 N'(k^2) ] + N(k^2)/( Sqrt[5] k^2 )
     */

     REAL parn = par->N->N(parktsqr, par->y);
     REAL n = par->N->N(ktsqr, par->y);

    if (std::abs(lnktsqr - par->lnktsqr) < 1e-15)
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
        result += ( n - std::exp(par->lnktsqr - lnktsqr)*parn )
            / std::abs( 1.0 - std::exp(par->lnktsqr - lnktsqr) );
    
    result += parn/std::sqrt( 4.0*std::exp(2.0*(lnktsqr - par->lnktsqr)) + 1.0);

    if (par->N->RunningCoupling() == MAXK)
        result *= Alphabar_s(std::max(ktsqr, parktsqr));
    else if (par->N->RunningCoupling() == MINK)
        result *= Alphabar_s(std::min(ktsqr, parktsqr));
    
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
REAL inthelperf_bkmom_constraint(REAL lnktsqr, void* p)
{
    ///TODO: Optimize to use lnktsqr and par->lnktsqr
    inthelper_bkmom* par = (inthelper_bkmom*) p;
    REAL ktsqr = std::exp(lnktsqr);
    REAL parktsqr = std::exp(par->lnktsqr);
    REAL n0 = par->N->N(parktsqr, par->y);  

    REAL result=0;
    if (std::abs(lnktsqr - par->lnktsqr) < 1e-15)
    {
        cerr << "ktsqr \\approx par->ktsqr and we can't handle this! y=" << par->y
            << " par->ktsqr=" << parktsqr << " at " << LINEINFO << endl;
            ///TODO
    }
    else
    {
        if (par->lnktsqr > lnktsqr)
            result += par->N->N(ktsqr, par->y);
        else
            result += par->N->N(ktsqr, par->y - (lnktsqr - par->lnktsqr));
        result -= std::exp(par->lnktsqr - lnktsqr)*n0;
        result /= std::abs( 1.0 - std::exp(par->lnktsqr - lnktsqr) );
    }
    result += n0/std::sqrt( 4.0*std::exp(2.0*(lnktsqr - par->lnktsqr)) + 1.0);

    if (par->N->RunningCoupling() == MAXK)
        result *= Alphabar_s(std::max(ktsqr, parktsqr));
    else if (par->N->RunningCoupling() == MINK)
        result *= Alphabar_s(std::min(ktsqr, parktsqr));

    return result;
}


REAL BruteForceSolver::RapidityDerivative(REAL ktsqr, REAL y)
{    
    inthelper_bkmom inthelp;
    inthelp.N=this;
    inthelp.y=y;
    inthelp.lnktsqr = std::log(ktsqr);
    gsl_function int_helper;

    if (kinematic_constraint==false)
        int_helper.function=&inthelperf_bkmom_noconstraint;
    else
        int_helper.function=&inthelperf_bkmom_constraint;
    int_helper.params=&inthelp;

    REAL result, abserr; 
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(KTSQRINTITERATIONS);

    REAL minlnktsqr = std::log(ktsqrvals[0]);
    REAL maxlnktsqr = std::log(ktsqrvals[ktsqrvals.size()-2]);

    int status;
    status=gsl_integration_qag(&int_helper, minlnktsqr,
            maxlnktsqr, 0, KTSQRINTACCURACY, 
            KTSQRINTITERATIONS, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    
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
 * TODO: Doesn't work, as RapidityDerivative doesn't use amplitud[] array
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
    // Find maxyind corresponding to maxy
    int maxyind=YPoints();
    for (unsigned int i=1; i<=YPoints(); i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }
    int largedifference=0;
    int starty=1;
    if (adams_method)
    {
        InitializeAdamsMethod();
        starty=2;
    }
    // ******************************
    for (int yind=starty; yind<=maxyind; yind++)
    {
        cout << "Solving for y=" << yvals[yind] << endl;
        // Solve N(y+DELTA_Y, kt) for every kt



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
            if (adams_method==true)
            {
                REAL old_der = derivatives[ktsqrind][yind-2];
                REAL adamsn = n[ktsqrind][yind-1] + dy*tmpder
                    + 1.0/2.0*dy*( tmpder - old_der);

                //cout << "reldiff at k=" << ktsqrvals[ktsqrind] <<": " << std::abs((newn-adamsn)/newn) << endl;
                newn = adamsn;
            }

               
                
            AddDataPoint(ktsqrind, yind, newn, tmpder );
            if( abs(newn - n[ktsqrind][yind-1])/n[ktsqrind][yind-1] > 0.1)
            {
                largedifference++;
            }

            /*if (ktsqrind % 100 == 0)
            {
                cout << "y=" << yvals[yind] << ", ktsqr=" << ktsqrvals[ktsqrind]
                    << " reldiff " << std::abs((newn-n[ktsqrind][yind-1])/n[ktsqrind][yind-1])
                    << " absdiff " << std::abs(newn - n[ktsqrind][yind-1]) << endl;
            }*/
            //cout << "N(ktsqr=" << ktsqrvals[ktsqrind] <<", y=" << yvals[yind] << ") = " << newn
            //<< ",  at lower rapidity it was " << n[ktsqrind][yind-1] << endl;
                
        }
    }
        

        
    cout << endl << "#" << largedifference << " out of " << maxyind * (KtsqrPoints()-1)
            << " too large differences" << endl;
    // Again
    if (averages>0) cerr << "Averagements are not well tested..." << endl;
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
    } 
    
    

}


/*
 * Intialize Adam's method
 * So we need to calculate amplitude at Y=yvals[1] using smaller step size,
 * then Solve() can use that amplitude and its derivative in 2nd order Adams
 * method
 */
void BruteForceSolver::InitializeAdamsMethod()
{
    n.clear();
    REAL y = yvals[1];
    yvals.clear();

    const int SMALLSTEPS = 10;

    for (uint i=0; i<=SMALLSTEPS; i++) // 10 steps to y=yvals[1]
    {
        yvals.push_back(y/static_cast<REAL>(SMALLSTEPS)*static_cast<REAL>(i));
    }

    for (uint i=0; i<=KtsqrPoints(); i++)   // Intialize every kt
    {
        std::vector<REAL> tmpvec;
        std::vector<REAL> tmpdervec;

        // y=0 initial condition
        tmpvec.push_back(InitialCondition( ktsqrvals[i] ));
        tmpdervec.push_back(0.0);
        for (unsigned int j=1; j<=YPoints(); j++)
        {
            tmpvec.push_back(0.0);
            tmpdervec.push_back(0.0);
        }
        n.push_back(tmpvec);
        derivatives.push_back(tmpdervec);
    }

    for (int yind=1; yind<=SMALLSTEPS; yind++)
    {
        cout << "Solving for y=" << yvals[yind] << endl;
        // Solve N(y+DELTA_Y, kt) for every kt
        #pragma omp parallel for
        for (int ktsqrind=0; ktsqrind<KtsqrPoints(); ktsqrind++)
        {
            REAL tmpkt = ktsqrvals[ktsqrind];
            REAL dy = yvals[yind]-yvals[yind-1];
            REAL tmpder = RapidityDerivative(tmpkt, yvals[yind-1]);
            REAL newn=n[ktsqrind][yind-1] + dy*tmpder;
            AddDataPoint(ktsqrind, yind, newn, tmpder );
        }
    }

    cout << " Amplitude at y=" << y << " solved using smaller step size, "
        << "starting to use larger step size and Adams' method." << endl;

    std::vector<REAL> amp; std::vector<REAL> der_y0, der_y1;
    for (uint i=0; i<=KtsqrPoints(); i++)
    {
        amp.push_back(n[i][SMALLSTEPS]);
        der_y0.push_back( derivatives[i][0] );
    }

    Initialize();

    for (uint i=0; i<=KtsqrPoints(); i++)
    {
        derivatives[i][0] = der_y0[i];
        n[i][1] = amp[i];
    }
    
}

BruteForceSolver::BruteForceSolver()
{
    rungekutta = false;
    adams_method=false;
}
BruteForceSolver::~BruteForceSolver()
{

}

void BruteForceSolver::SetAdamsMethod(bool ad)
{
    adams_method = ad;
}

void BruteForceSolver::SetRungeKutta(bool rk)
{
    rungekutta = rk;
}
