/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
#include <tools/interpolation.hpp>
#include "solver_force.hpp"
#include <algorithm>
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
    REAL lnktsqr;
    BruteForceSolver* N;
    REAL y;
    REAL offset;
    const REAL* array;    // Array containing values used in Runge Kutta method
    bool rungekutta;
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

    REAL parn, n;

    if (par->rungekutta)
    {
        parn = par->N->InterpolateN(par->lnktsqr, par->array);
        n = par->N->InterpolateN(lnktsqr, par->array);
    } else
    {
        parn = par->N->N(parktsqr, par->y);
        n = par->N->N(ktsqr, par->y);
    }

    if (std::abs(lnktsqr - par->lnktsqr) < 1e-5)
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
        result += (n - parktsqr/ktsqr * parn) / std::abs( 1.0 - parktsqr/ktsqr);
      
    //    result += ( n - std::exp(par->lnktsqr - lnktsqr)*parn )
    //        / std::abs( 1.0 - std::exp(par->lnktsqr - lnktsqr) );

    //result += parn/std::sqrt( 4.0*std::exp(2.0*(lnktsqr - par->lnktsqr)) + 1.0);
    result += parn/std::sqrt( 4.0*SQR(ktsqr / parktsqr) + 1.0 );

    if (par->N->RunningCoupling() == MAXK)
        result *= Alphabar_s(std::max(ktsqr, parktsqr), par->N->AlphasScaling());
    else if (par->N->RunningCoupling() == MINK)
        result *= Alphabar_s(std::min(ktsqr, parktsqr), par->N->AlphasScaling());

    if (isnan(result))
        cerr << "error at ktsqr=" << ktsqr << " parktsqr " << parktsqr << " n " << n
        << " parn " << parn << " res " << result << endl;
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
    if (std::abs(lnktsqr - par->lnktsqr) < 1e-5)
    {
        //cerr << "ktsqr \\approx par->ktsqr and we can't handle this! y=" << par->y
        //    << " par->ktsqr=" << parktsqr << " at " << LINEINFO << endl;
            ///TODO
    }
    else
    {
        if (par->lnktsqr > lnktsqr)
            result += par->N->N(ktsqr, par->y);
        else
            result += par->N->N(ktsqr, par->y - (lnktsqr - par->lnktsqr));
        //result -= std::exp(par->lnktsqr - lnktsqr)*n0;
        result -= parktsqr / ktsqr * n0;
        //result /= std::abs( 1.0 - std::exp(par->lnktsqr - lnktsqr) );
        result /= std::abs( 1.0 - parktsqr/ktsqr );
    }
    //result += n0/std::sqrt( 4.0*std::exp(2.0*(lnktsqr - par->lnktsqr)) + 1.0);
    result += n0/std::sqrt( 4.0* SQR(ktsqr/parktsqr) + 1.0 );

    if (par->N->RunningCoupling() == MAXK)
        result *= Alphabar_s(std::max(ktsqr, parktsqr), par->N->AlphasScaling());
    else if (par->N->RunningCoupling() == MINK)
        result *= Alphabar_s(std::min(ktsqr, parktsqr), par->N->AlphasScaling());

    if (isnan(result) or isinf(result))
    {
        cerr << "Integrand is " << result << " at " << LINEINFO <<", ktsqr=" <<ktsqr <<
        " n0 " << n0 << " parktsqr " << parktsqr << " y " << par->y << " real_n ";
        if (par->lnktsqr > lnktsqr)
            cerr << "local " << par->N->N(ktsqr, par->y);
        else
             cerr << "nonlocal " << par->N->N(ktsqr, par->y - (lnktsqr - par->lnktsqr))
            << " dy " << lnktsqr-par->lnktsqr;
        cerr << endl;
        result=0;
        exit(1);
    }
    
    return result;
}

/*
 * Calculate \partial_Y N
 * If array != NULL (default), then we are using runge kutta and this array
 * should be used as a N(ktsqr) (doesn't work with kinematicla constraint)
 */

REAL BruteForceSolver::RapidityDerivative(REAL ktsqr, REAL y, const REAL* array)
{
    inthelper_bkmom inthelp;
    inthelp.N=this;
    inthelp.y=y;
    inthelp.lnktsqr = std::log(ktsqr);
    inthelp.array = array;
    inthelp.rungekutta=rungekutta;
    gsl_function int_helper;

    if (kinematic_constraint==false)
        int_helper.function=&inthelperf_bkmom_noconstraint;
    else
        int_helper.function=&inthelperf_bkmom_constraint;
    int_helper.params=&inthelp;

    REAL result, abserr;
    //int ktsqriter = KTSQRINTITERATIONS;
    int ktsqriter = KtsqrPoints();
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(ktsqriter);

    REAL minlnktsqr = std::log(ktsqrvals[0]);
    REAL maxlnktsqr = std::log( Ktsqrval(KtsqrPoints()-1) );

    int status;
    status=gsl_integration_qag(&int_helper, minlnktsqr,
            maxlnktsqr, 0, KTSQRINTACCURACY, 
            ktsqriter, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
    
    gsl_integration_workspace_free(workspace);
    if (status) { cerr << "Error " << status << " at " << LINEINFO << ":"
        << " ktsqr=" << ktsqr <<", y=" << y << " result=" << result << ", abserror=" <<
        abserr << " relerror: " << abserr/result << endl; }

    // Nonlinear term
    // if rc is {MIN,MAX}k, then result is already multiplied by alpha_s
    
    if (RunningCoupling()==MAXK or RunningCoupling() == MINK)
        result -= SQR(N(ktsqr, y))*Alphabar_s(ktsqr, alphas_scaling);
    else
        result -= SQR(N(ktsqr, y));

    if (RunningCoupling()==CONSTANT)
        result*=ALPHABAR_s;
    if (RunningCoupling()==PARENT_DIPOLE)
        result *= Alphabar_s(ktsqr, alphas_scaling);
    
    return result;
}
/*
 * GSL ODEIV evolution
 * Returns RapidityDerivative for every ktsqr in vector
 *
 */
struct EvolutionHelper
{
    BruteForceSolver* N;
};
int Evolve(REAL y, const REAL amplitude[], REAL dydt[], void *params)
{
    EvolutionHelper* helper = (EvolutionHelper*)params;
   #pragma omp parallel for
    for( int i=0 ; i<helper->N->KtsqrPoints() ; i++) {
        REAL res = helper->N->RapidityDerivative(helper->N->Ktsqrval(i), y, amplitude);
        dydt[i] = res;
   }
   return GSL_SUCCESS;
}

const REAL MINSTEP=0.002;
//const REAL MINSTEP=0.001;
// Solve BK, lowest order
void BruteForceSolver::Solve(REAL maxy)
{

    REAL Y=0;   // We are allways solved amplitude up to Y
    int yind=0; // yvals[yind]=Y always
    REAL nexty=delta_y/4;   // smaller step size a the beginning

    // **** used in GSL solver *****
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;

    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, KtsqrPoints());
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0, 0.01);    //abserr relerr
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (KtsqrPoints());
    EvolutionHelper help; help.N=this;
    gsl_odeiv_system sys = {Evolve, NULL, KtsqrPoints(), &help};
    REAL Yi=0;
    REAL h = delta_y;   // Original step size
    ///////

    // these arrays are used to temporarily save solved amplitudes
    REAL *amplitude=new REAL[KtsqrPoints()];
    for (int ktsqrind=0; ktsqrind < KtsqrPoints(); ktsqrind++)
    {
        amplitude[ktsqrind] = n[0][ktsqrind] ;  
    }

    REAL *ders = new REAL[KtsqrPoints()];

    
    if (adams_method)
    {
        cerr << "Adams method is not working with new system... FIXME!" << endl;
        exit(1);
        //InitializeAdamsMethod();
        //starty=2;
    }
    // ******************************
    do
    {
        
        // Solve N(y+DELTA_Y, kt) for every kt

        if (rungekutta)
        {
            cout << "Solving for y=" << nexty << endl;
            Yi = nexty;
            // gsl_odeiv_evolve_apply increases Y according to the step size
            while (Y<Yi)
            {
                int status = gsl_odeiv_evolve_apply(e, c, s, &sys,
                        &Y, Yi, &h, amplitude);
                if (status != GSL_SUCCESS) {
                    cerr << "Error in gsl_odeiv_evolve_apply at " << LINEINFO
                        << ": " << gsl_strerror(status) << " (" << status << ")"
                        << " y=" << Y << ", h=" << h << endl;
                }
            } 

            AddRapidity(nexty); // During the evolution Y has evolved up to nexty  
            for (int ktsqrind = 0; ktsqrind < KtsqrPoints(); ktsqrind++)
            {
                if (ktsqrind+1< KtsqrPoints())
                {
                    if (amplitude[ktsqrind+1] > amplitude[ktsqrind])
                    cerr << "Amplitude increases as k increases! index " << ktsqrind+1 << "/" << KtsqrPoints()
                    << " " << amplitude[ktsqrind] << "->" << amplitude[ktsqrind+1] << endl;
                } 
                AddDataPoint(ktsqrind, yind+1, amplitude[ktsqrind], 0.0);
            }

            nexty += delta_y;
            yind++;
            continue;
        }
        // If we end up here, we are not using Runge Kutta
        REAL dy = nexty-yvals[yind];
        bool ok=true;
        #pragma omp parallel for
        for (int ktsqrind=0; ktsqrind<KtsqrPoints(); ktsqrind++)
        {
            if (!ok) continue;  // We will decrease the step size
            REAL tmpktsqr = ktsqrvals[ktsqrind];
            REAL tmpder = RapidityDerivative(tmpktsqr, yvals[yind]);
            REAL oldn = n[yind][ktsqrind];
            REAL newn= oldn + dy*tmpder;

            // Adaptive step size (quite stupid one)
            // If reldiff > 0.1 and absdiff > 1e-5 then use smaller step size
            //if (std::abs(dy*tmpder) > 1e-3 and std::abs(dy*tmpder/oldn > 0.01)
            //    and (nexty-Y)*2.0/3.0 > MINSTEP)
            
            #pragma omp critical
            {
                if (ok==true and std::abs(dy*tmpder/oldn)>0.01 and (nexty-Y)*2.0/3.0 > MINSTEP)
                {
                    nexty = Y + (nexty-Y)*2.0/3.0;
                    cout << "Taking smaller step size " <<nexty-Y <<". ";
                    ok=false;
                }
                else
                {
                    amplitude[ktsqrind]=newn;
                    ders[ktsqrind] = tmpder;
                }
            }
        }

        if (!ok)
            continue;

        // Ok, solved for all k_T withour errors, save
        AddRapidity(nexty);
        REAL step = nexty-Y;
        for (uint i=0; i<KtsqrPoints(); i++)
        {
            if (i+1<KtsqrPoints())
            {
                if (amplitude[i+1] > amplitude[i])
                {
                    cerr << "Amplitude increases as k increases! index " << i+1 << "/" << KtsqrPoints()
                    << " " << amplitude[i] << "->" << amplitude[i+1] << endl;
                }
            }
            AddDataPoint(i, yind+1, amplitude[i], ders[i] );
        }
        cout << "Solved y=" << nexty << " step " << step << endl;
        
        if (1.5*step < delta_y)
            step*=1.5;
        else
           step = delta_y;
        Y = nexty;
        nexty += step;
        yind++;
    }while(Y <= maxy+delta_y);

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    delete[] amplitude;
    delete[] ders;
}


/*
 * Intialize Adam's method
 * So we need to calculate amplitude at Y=yvals[1] using smaller step size,
 * then Solve() can use that amplitude and its derivative in 2nd order Adams
 * method
 */
void BruteForceSolver::InitializeAdamsMethod()
{
    cerr << "Adam's method not implemented!" << endl;
}

/*
 * InterpolateN
 * Calculates amplitude at given ktsqr using the given array
 * array[i] = amplitude at ktsqr=ktsqrvals[i]
 * If bspline is true, use bspline, otherwise spline
 * If der is true, return derivative (TODO: not implemented)
 * By default bspline and der = false
 */
REAL BruteForceSolver::InterpolateN(REAL lnktsqr, const REAL* array, bool bspline, bool der)
{
    if (lnktsqr >= lnktsqrvals[KtsqrPoints()-1]) return array[KtsqrPoints()-1];
    if (lnktsqr <= lnktsqrvals[0]) return array[0];
    REAL kt = std::exp(0.5*lnktsqr);
        
    int ktsqrind = FindIndex(lnktsqr, lnktsqrvals);
    if (ktsqrind < 0)
        cerr << "Negative ktsqrindex at " << LINEINFO << endl;


    // Keep y fixed, interpolate ktsqr
    // Interpolate only INTERPOLATION_POINTS points in order to make this
    // almost fast

    unsigned int interpolation_start, interpolation_end;
    
    if (ktsqrind - interpolation_points/2 < 0)
    {
		interpolation_start=0;
		interpolation_end=interpolation_points;
	}
	else if (ktsqrind + interpolation_points/2 > KtsqrPoints()-1 )
	{
		interpolation_end = KtsqrPoints()-1;
		interpolation_start = KtsqrPoints()-interpolation_points-2;
	}
	else
	{
		interpolation_start = ktsqrind - interpolation_points/2;
		interpolation_end = ktsqrind + interpolation_points/2;
	}
	int interpo_points = interpolation_end - interpolation_start+1;
    
    REAL *tmparray = new REAL[interpo_points];
    REAL *tmpxarray = new REAL[interpo_points];
    for (int i=interpolation_start; i<= interpolation_end; i++)
    {
		tmpxarray[i-interpolation_start]=ktvals[i];
        tmparray[i-interpolation_start] = array[i];	
    }

    Interpolator interp(tmpxarray, tmparray, interpo_points);
    if (bspline)
        interp.SetMethod(INTERPOLATE_BSPLINE);
    interp.Initialize();
    REAL res = interp.Evaluate(kt);

    delete[] tmparray;
    delete[] tmpxarray;
    return res;
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
