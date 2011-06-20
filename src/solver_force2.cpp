/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
#include "solver_force2.hpp"
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <cmath>
using std::abs;

const REAL UINTACCURACY = 0.001;
const int UINTITERATIONS = 10000;

// Amplitude class methods which solve BK directly (=brute force)


/*
 * \partial_Y from BK in momentum space
 */

struct inthelper_bkmom2
{
    REAL u;
    BruteForceSolver2* N;
    REAL y;
};


/*
 * BK equation in momentum space integrated over \theta (b-indep. situation)
 * Ref e.g. hep-ph/0110325 eq (8)
 */
REAL inthelperf_bkmom2_noconstraint(REAL v, void* p)
{
    inthelper_bkmom2* par = (inthelper_bkmom2*) p;

    REAL result=0;

    REAL nu = par->N->UAmplitude(par->u, par->y);

    /* When v -> par->u first doesn't diverge, but numerics would fail
     */

    if (v > 1.0-1e-45 or par->u>1.0-1e-45) return 0;

    if (std::abs(par->u - v) < 1e-45)
    {
        // Let's say that the first term is zero
    }
    else   
        result = (v/(1.0-v)*par->N->UAmplitude(v, par->y) - par->u/(1.0-par->u)*nu) /
            std::abs( v/(1.0-v) - par->u/(1.0-par->u) );
    

    result += par->u/(1.0-par->u) * nu
        / std::sqrt(4.0*SQR( v/(1.0-v) ) + SQR( par->u / (1.0-par->u) ) );

    result /= v*(1.0-v);
    
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
REAL inthelperf_bkmom2_constraint(REAL ktsqr, void* p)
{
    /* inthelper_bkmom* par = (inthelper_bkmom*) p;

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

    return result; */
    return 0.0;
}

REAL BruteForceSolver2::RapidityDerivative(REAL u, REAL y)
{
    REAL alphabar=0.2;
    if (RunningCoupling())
        alphabar = Alpha_s(Ktsqr(u))*Nc/M_PI;
    
    inthelper_bkmom2 inthelp;
    inthelp.N=this;
    inthelp.y=y;
    inthelp.u = u;
    gsl_function int_helper;

    if (kinematic_constraint==false)
        int_helper.function=&inthelperf_bkmom2_noconstraint;
    else
        int_helper.function=&inthelperf_bkmom2_constraint;
    int_helper.params=&inthelp;

    REAL result, abserr; 
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(UINTITERATIONS);

    int status;
    // If ktsqr != ktsqr', we can tell GSL that there is a difficult point
    // at ktsqr'=ktsqr
   /* if (ktsqr>ktsqrvals[0] and ktsqr<ktsqrvals[ktsqrvals.size()-2] and true==false)
    {
        REAL range[3]; range[0]=ktsqrvals[0];
        range[1]=ktsqr; range[2]=ktsqrvals[ktsqrvals.size()-1];
        status = gsl_integration_qagp(&int_helper, range, 3, 0, KTSQRINTACCURACY,
            KTSQRINTITERATIONS, workspace, &result, &abserr);
    } else*/
    {
    status=gsl_integration_qag(&int_helper, minu, 1.0, 0, UINTACCURACY, 
             UINTITERATIONS, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    }
    gsl_integration_workspace_free(workspace);
    if (status ) cerr << "Error " << status << " at " << LINEINFO << ":"
        << " u=" << u <<", y=" << y << " result=" << result << ", abserror=" <<
        abserr << " relerror: " << abserr/result << endl;

    // Nonlinear term
    result -= SQR(UAmplitude(u, y) );
    
    return alphabar*result;
}


// Solve BK, lowest order
void BruteForceSolver2::Solve(REAL maxy)
{
    Prepare();
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
        for (int uind=0; uind<upoints-1; uind++)
        {
            REAL tmpu = uvals[uind];
            REAL dy = yvals[yind]-yvals[yind-1];
            REAL tmpder = RapidityDerivative(tmpu, yvals[yind-1]);


            REAL newn=amplitude[uind][yind-1] + dy*tmpder;

            // Adams method: apprximate derivative as a 2nd order polynomial
            // Y_{i+1} = y_i + hf_i + 1/2 h (f_i - f_{i-1} )
            // We can use this only for yind>1
            /*if (yind>1 and adams_method==true)
            {
                REAL old_der = derivatives[ktsqrind][yind-2];
                newn = amplitude[uind][yind-1] + dy*tmpder
                    + 1.0729400/2.0*dy*( tmpder - old_der);
            }
            */
            
            
            amplitude[uind][yind] = newn;
            ktsqrvals[uind] = Ktsqr(uvals[uind]);
            n[uind][yind] = newn;
            if( abs(newn - amplitude[uind][yind-1])/amplitude[uind][yind-1] > 0.1)
            {
                largedifference++;
            }
            //cout << "N(ktsqr=" << ktsqrvals[ktsqrind] <<", y=" << yvals[yind] << ") = " << newn
            //<< ",  at lower rapidity it was " << n[ktsqrind][yind-1] << endl;
            
        }
    }
    cout << endl << "#" << largedifference << " out of " << maxyind * (KtsqrPoints()-1)
            << " too large differences" << endl;

}

/*
 * Amplitude as a function of u and y
 * Sets explicitly N(u)=0
 */
REAL BruteForceSolver2::UAmplitude(REAL u, REAL y)
{
    if (y<eps and datafile==false) return InitialCondition(Ktsqr(u));
    if (y<eps) y=0;
    
    // Find uval and yval indexes refer to index for which
    // val[index]  is smaller than (or equal) y/ktsqr
    int yind = -1;
    int uind=-1;
    for (unsigned int i=0; i<yvals.size()-1; i++)
    {
        if (yvals[i]<=y and yvals[i+1]>y)
        {
            yind=i;
            break;
        }
    }
    if (yind < 0) // Didn't find, so refers to the largest one 
    {
        yind=yvals.size()-1;
        if (y-yvals[yind]>0.05)
            cerr << "Asked amplitude at too large Y=" << y << ", falling back to "
                << " y=" << yvals[yind] << ". " << LINEINFO << endl;
    }
    for (unsigned int i=0; i<uvals.size()-1; i++)
    {
        if (uvals[i]<=u and uvals[i+1]>u)
        {
            uind=i;
            break;
        }
    }
    if (uind < 0) // Didn't, shouldn't end up here
    {
        cerr << "Didn't find uind for u=" << u << " at " << LINEINFO << endl;
    }

    // Keep y fixed, interpolate ktsqr
    // Interpolate only INTERPOLATION_POINTS points in order to make this
    // almost fast
    // Interpolate linearly in y and use spline in ktsqr

    unsigned int interpolation_start, interpolation_end;
    if (uind - INTERPOLATION_POINTS/2 < 0)
    {
		interpolation_start=0;
		interpolation_end=INTERPOLATION_POINTS;
	}
	else if (uind + INTERPOLATION_POINTS/2 > upoints-1 )
	{
		interpolation_end = upoints-1;
		interpolation_start = upoints-INTERPOLATION_POINTS-1;
	}
	else
	{
		interpolation_start = uind - INTERPOLATION_POINTS/2;
		interpolation_end = uind + INTERPOLATION_POINTS/2;
	}
	int interpo_points = interpolation_end - interpolation_start+1;
    
    REAL *tmparray = new REAL[interpo_points];
    REAL *tmpxarray = new REAL[interpo_points];
    for (uint i=interpolation_start; i<= interpolation_end; i++)
    {
		tmpxarray[i-interpolation_start]=uvals[i];

        tmparray[i-interpolation_start] = amplitude[i][yind];

        // Interpolate in y if possible
		if (yind < yvals.size()-1 )
        {
            if (amplitude[i][yind+1]>-1e-20 and y-yvals[yind]>0.000001 )
            {
                tmparray[i-interpolation_start]=amplitude[i][yind] 
                 + (y - yvals[yind])* (amplitude[i][yind+1] - amplitude[i][yind])
                    / (yvals[yind+1]-yvals[yind]);
            }
		} 
		
			
    }
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, interpo_points);
    gsl_spline_init(spline, tmpxarray, tmparray, interpo_points);

    REAL res;
    int status = gsl_spline_eval_e(spline, u, acc, &res);
    if (status)
    {
        cerr << "Interpolatioin failed at " << LINEINFO << ", error " << gsl_strerror(status)
         << " (" << status << "), u=" << u << ", y=" << y << endl;
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    delete[] tmparray;
    delete[] tmpxarray;

    return res;
    

}

REAL BruteForceSolver2::Ktsqr(REAL u)
{
    return u/(1-u);
}

REAL BruteForceSolver2::U(REAL ktsqr)
{
    return ktsqr/(1.0+ktsqr);
}

BruteForceSolver2::BruteForceSolver2()
{

}

void BruteForceSolver2::Prepare()
{
    // Intialize n[u][y]
    upoints = KtsqrPoints();
    minu = U(minktsqr);
    maxu = U(maxktsqr);

    for (uint i=0; i<=upoints-2; i++)
    {
        uvals.push_back(minu + (maxu-minu)/(upoints-1)*static_cast<REAL>(i) );
    }
    uvals.push_back(1.0);
    
    for (unsigned int i=0; i<=upoints-2; i++)   // Intialize every u
    {
        std::vector<REAL> tmpvec;

        // y=0 initial condition
        tmpvec.push_back(InitialCondition( Ktsqr(uvals[i]) ) );
        for (unsigned int j=1; j<=YPoints(); j++)
        {
            tmpvec.push_back(0.0);
        }
        amplitude.push_back(tmpvec);
    }

    // Last case, u=1, N=0
    std::vector<REAL> tmpvec;
    for (uint j=0; j<=YPoints(); j++)
        tmpvec.push_back(0.0);
    amplitude.push_back(tmpvec);
}

BruteForceSolver2::~BruteForceSolver2()
{

}
