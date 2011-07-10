/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
#include "interpolation.hpp"
#include "solver_force.hpp"
#include <algorithm>
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
        parn = par->N->InterpolateN(parktsqr, par->array);
        n = par->N->InterpolateN(ktsqr, par->array);
    } else
    {
        parn = par->N->N(parktsqr, par->y);
        n = par->N->N(ktsqr, par->y);
    }

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
    int ktsqriter = KTSQRINTITERATIONS;
    if (ktsqriter >= KtsqrPoints()) ktsqriter = KtsqrPoints()-1;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(ktsqriter);

    REAL minlnktsqr = std::log(ktsqrvals[0]);
    REAL maxlnktsqr = std::log( Ktsqrval(KtsqrPoints()-1) );

    int status;
    status=gsl_integration_qag(&int_helper, minlnktsqr,
            maxlnktsqr, 0, KTSQRINTACCURACY, 
            ktsqriter, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    
    gsl_integration_workspace_free(workspace);
    if (status ) cerr << "Error " << status << " at " << LINEINFO << ":"
        << " ktsqr=" << ktsqr <<", y=" << y << " result=" << result << ", abserror=" <<
        abserr << " relerror: " << abserr/result << endl;

    // Nonlinear term
    // if rc is {MIN,MAX}k, then result is already multiplied by alpha_s
    if (RunningCoupling()==MAXK or RunningCoupling() == MINK)
        result -= SQR(N(ktsqr, y))*Alphabar_s(ktsqr);
    else
        result -= SQR(N(ktsqr, y));

    if (RunningCoupling()==CONSTANT)
        result*=ALPHAS;
    if (RunningCoupling()==PARENT_DIPOLE)
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
int Evolve(REAL y, const REAL amplitude[], REAL dydt[], void *params)
{
    cout << "Evolving with y=" << y << endl;
    EvolutionHelper* helper = (EvolutionHelper*)params;
   // Actually we don't need amplitude[] as it only contains
   // elements N[ktsqrind][yind], but GSL uses it to determine whether
   // the volution is accurate 
   #pragma omp parallel for
   for( int i=0 ; i<helper->N->KtsqrPoints() ; i++) {
      dydt[i] = helper->N->RapidityDerivative(helper->N->Ktsqrval(i), y, amplitude);
   }
   return GSL_SUCCESS;
}


// Solve BK, lowest order
void BruteForceSolver::Solve(REAL maxy)
{

    int largedifference=0;
    int starty=1;
    REAL Y=0;   // We are allways solved amplitude up to Y
    int yind=0; // yvals[yind]=Y always
    REAL nexty=delta_y;

    // **** used in GSL solver *****
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;

    gsl_odeiv2_step * s    = gsl_odeiv2_step_alloc (T, KtsqrPoints());
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (0.0, 0.05);    //abserr relerr
    gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (KtsqrPoints());
    EvolutionHelper help; help.N=this;
    gsl_odeiv2_system sys = {Evolve, NULL, KtsqrPoints(), &help};
    REAL Yi=0;
    REAL h = delta_y;   // Original step size
    ///////

    // these arrays are used to temporarily save solved amplitudes
    REAL *amplitude=new REAL[KtsqrPoints()];
    for (int ktsqrind=0; ktsqrind < KtsqrPoints(); ktsqrind++)
    {
        amplitude[ktsqrind] = std::exp( ln_n[0][ktsqrind] );   // TODO: log?
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
        cout << "Solving for y=" << nexty << endl;
        // Solve N(y+DELTA_Y, kt) for every kt

        if (rungekutta)
        {
            Yi = nexty;
            // gsl_odeiv_evolve_apply increases Y according to the step size
            // why gsl_odeiv_evolve_apply doesn't set Y=Y_i at the end?
            while (Y<Yi)
            {
                int status = gsl_odeiv2_evolve_apply(e, c, s, &sys,
                        &Y, Yi, &h, amplitude);
                if (status != GSL_SUCCESS) {
                    cerr << "Error in gsl_odeiv_evolve_apply at " << LINEINFO
                        << ": " << gsl_strerror(status) << " (" << status << ")"
                        << " y=" << Y << ", h=" << h << endl;
                }
                cout << "Evolved up to " << Y << "/" << Yi << ", h=" << h << endl;
            } // end while (useless loop?)
            cout << "Solved yind " << yind << " to Y=" << Y << " with step size "
                    << h << endl;
            AddRapidity(nexty); // During the evolution Y has evolved up to nexty  
            for (int ktsqrind = 0; ktsqrind < KtsqrPoints(); ktsqrind++)
            {
                AddDataPoint(ktsqrind, yind+1, amplitude[ktsqrind], 0.0);
                //n[ktsqrind][yind] = amplitude[ktsqrind];
            }

            nexty += delta_y;
            yind++;
            continue;
        }
        
        // If we end up here, we are not using Runge Kutta

        #pragma omp parallel for
        for (int ktsqrind=0; ktsqrind<KtsqrPoints(); ktsqrind++)
        {
            REAL tmpkt = ktsqrvals[ktsqrind];
            REAL dy = nexty-yvals[yind];
            REAL tmpder = RapidityDerivative(tmpkt, yvals[yind]);
            REAL newn=std::exp(ln_n[yind][ktsqrind]) + dy*tmpder;

            // Adams method: apprximate derivative as a 2nd order polynomial
            // Y_{i+1} = y_i + hf_i + 1/2 h (f_i - f_{i-1} )
            // We can use this only for yind>1
            if (adams_method==true)
            {
                REAL old_der = derivatives[yind-1][ktsqrind];
                REAL adamsn = std::exp(ln_n[yind][ktsqrind] ) + dy*tmpder
                    + 1.0/2.0*dy*( tmpder - old_der);

                //cout << "reldiff at k=" << ktsqrvals[ktsqrind] <<": " << std::abs((newn-adamsn)/newn) << endl;
                newn = adamsn;
            }

            if( abs(newn - std::exp(ln_n[yind][ktsqrind]) ) / std::exp(ln_n[yind][ktsqrind]) > 0.2)
            {
                largedifference++;
            }

            amplitude[ktsqrind]=newn;
            ders[ktsqrind] = tmpder;
        }

        // Ok, solved for all k_T withour errors, save
        AddRapidity(nexty);
        for (uint i=0; i<KtsqrPoints(); i++)
        {
            AddDataPoint(i, yind+1, amplitude[i], ders[i] );
        }
        Y = nexty;
        nexty += delta_y;
        yind++;
    }while(Y < maxy);
        

        
    cout << endl << "#" << largedifference << " too large differences" << endl;
    // Again
    /*
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
    */

    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
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
    cerr << "TODO at " << LINEINFO << ": adam's method is not impelmented with "
    << "tabulated ln_n" << endl;
    /*
    n.clear();
    REAL y = yvals[1];
    yvals.clear();
    derivatives.clear();

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
    */
    
}

/*
 * InterpolateN
 * Calculates amplitude at given ktsqr using the given array
 * array[i] = amplitude at ktsqr=ktsqrvals[i]
 * If bsplien is true, use bspline, otherwise spline
 * If der is true, return derivative (TODO: not implemented)
 * By default bspline and der = false
 */
REAL BruteForceSolver::InterpolateN(REAL ktsqr, const REAL* array, bool bspline, bool der)
{
    int ktsqrind = FindIndex(ktsqr, ktsqrvals);
    if (ktsqr >= MaxKtsqr()) // Didn't find, so refers to the largest one 
        return array[ktsqrvals.size()-1];
    if (ktsqr <= MinKtsqr())
        return array[0];


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
		interpolation_end = KtsqrPoints();
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
		tmpxarray[i-interpolation_start]=ktsqrvals[i];
        tmparray[i-interpolation_start] = array[i];	
    }

    Interpolator interp(tmpxarray, tmparray, interpo_points);
    if (bspline)
        interp.SetMethod(INTERPOLATE_BSPLINE);
    interp.Initialize();
    REAL res = interp.Evaluate(ktsqr);

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
