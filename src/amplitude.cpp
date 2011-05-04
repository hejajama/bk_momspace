/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
//#include <bci.h>
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <cmath>


using std::abs;

Amplitude::Amplitude()
{
    for (int i=0; i<POINTS_KTSQR; i++)
        ktsqrvals.push_back(MINKTSQR * std::pow(KTSQR_MULTIPLIER, i) );
    for (int i=0; i<POINTS_Y; i++)
        yvals.push_back((REAL) i * DELTA_Y);

}

/*
 * Initialize
 * Initial condition must be set
 */
void Amplitude::Initialize()
{
    for (int i=0; i<POINTS_KTSQR; i++)   // Intialize every kt
    {
        std::vector<REAL> tmpvec;
        std::vector<REAL> tmpdervec;

        // y=0 initial condition
        tmpvec.push_back(InitialCondition( ktsqrvals[i] ));
        tmpdervec.push_back(-1.0);
        for (int j=1; j<POINTS_Y; j++)
        {
            tmpvec.push_back(-1.0);
            tmpdervec.push_back(-1.0);
        }
        n.push_back(tmpvec);
        derivatives.push_back(tmpdervec);
    }
}

/*
 * Compute amplitude from the tabulated values
 */

REAL Amplitude::N(REAL ktsqr, REAL y)
{
    if (y<eps) return InitialCondition(ktsqr);

    // Linear interpolation or even extrapolation
    // TODO: SPLINE

    // Find ktsqrval and yval indexes refer to index for which
    // val[index]  is smaller than (or equal) y/ktsqr
    int yind = -1;
    int ktsqrind=-1;
    for (int i=0; i<=yvals.size()-2; i++)
    {
        if (yvals[i]<=y and yvals[i+1]>y)
        {
            yind=i;
            break;
        }
    }
    if (yind < 0) // Didn't find, so refers to the largest one 
    {
        yind=yvals.size()-2;
        cerr << "Asked amplitude at too large Y, falling back to "
            << " y=" << yvals[yind] << ". " << LINEINFO << endl;
    }
    for (int i=0; i<=ktsqrvals.size()-2; i++)
    {
        if (ktsqrvals[i]<=ktsqr and ktsqrvals[i+1]>ktsqr)
        {
            ktsqrind=i;
            break;
        }
    }
    if (ktsqrind < 0) // Didn't find, so refers to the largest one 
    {
        ktsqrind=ktsqrvals.size()-2;
        cerr << "Asked amplitude at too large ktsqr=" << ktsqr << ", falling back to "
            << " ktsqr=" << ktsqrvals[ktsqrind] << ". " << LINEINFO << endl;
    }

    // Keep y fixed, interpolate ktsqr
    ///FIXME: SSLLLLOOOOWW
    /*
    REAL *tmparray = new double[POINTS_KTSQR];
    REAL *tmpxarray = new double[POINTS_KTSQR];
    for (int i=0; i< POINTS_KTSQR; i++)
    {
        tmparray[i]=n[i][yind];
        tmpxarray[i]=ktsqrvals[i];
    }
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, POINTS_KTSQR);
    gsl_spline_init(spline, tmpxarray, tmparray, POINTS_KTSQR);

    REAL res = gsl_spline_eval(spline, ktsqr, acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    ///TODO: interpolate in y
    */
    // Linear
    REAL linear = n[ktsqrind][yind] +
        (ktsqr - ktsqrvals[ktsqrind]) *
            (n[ktsqrind+1][yind] - n[ktsqrind][yind]) / (ktsqrvals[ktsqrind+1] - ktsqrvals[ktsqrind]);
    if (n[ktsqrind][yind+1]>-eps)   // We can also interpolate y
        linear += (y - yvals[yind]) * (n[ktsqrind][yind+1] - n[ktsqrind][yind]) / (yvals[yind+1]-yvals[yind]);
    

    //cout << "Spline: " << res << "  linear: " << linear << endl;
    return linear;

}

/*
 * Interpolate
 * Interpolates data using libbci (cubic 2d interpolation)
 * Assumes that n[ktsqr][y] is filled up to y=y_max
 * Performance probably wery bad
 *
 * TODO: REQUIRES SQUARE DATA MATRIX :(
 */
void Amplitude::Interpolate()
{/*
    if (interpolated_amplitude != NULL)
        delete[] interpolated_amplitude;

    // Search largest yind
    int maxyind=-1;
    for (int i=POINTS_Y-1; i>=0; i--)
    {
        if (n[0][i]>-eps)
        {
            maxyind=i;
            break;
        }
    }
    if (maxyind==-1)
    {
        cerr << "Didn't find maxyind! " << LINEINFO << endl;
        return;
    }

    cout << "Interpolating up to yind=" << maxyind << endl;

    doublexyz *tmparray = new doublexyz[POINTS_KTSQR * (maxyind+1)];
    for (int ktind = 0; ktind < POINTS_KTSQR; ktind++)
    {
        for (int yind=0; yind <= maxyind; yind++)
        {
            tmparray[maxyind*ktind+yind].x = ktsqrvals[ktind];
            tmparray[maxyind*ktind+yind].y = yvals[yind];
            tmparray[maxyind*ktind+yind].z = n[ktind][yind];
        }
    }

    // 10 times more points in x and y axis
    interpolated_amplitude = new doublexyz[POINTS_KTSQR*10*(maxyind+1)*10];
    td_fillgrid(tmparray, POINTS_KTSQR, maxyind+1, interpolated_amplitude,
        POINTS_KTSQR*10, (maxyind+1)*10 );

    delete[] tmparray;
*/    
    

}

/*
 * Add computed data point, e.g. N(ktsqr, y), where ktsqr=ktsqrvals[ktsqrindex]
 * etc. der is \partial_Y N( ktsqr, yvals[yindex-1] ), e.g. the rapidity
 * derivative used to compute value
 */

void Amplitude::AddDataPoint(int ktsqrindex, int yindex, REAL value, REAL der)
{
    // NB: Assumes that we allways move higher in rapidity/ktsqr: when we solve
    // BK at Y=Y_0+DELTA_Y, we start at ktsqrindex=0 and move to ktsqrindex=max
    // keeping Y=Y_0+DELTA_Y, then start at ktsqrindex=0 and Y=Y_0+2DELTA_Y

    if (ktsqrindex < 0 or yindex < 0 or ktsqrindex > n.size()-1
        or yindex > n[ktsqrindex].size()-1)
    {
        cerr << "Invalid ktsqr/y index " << ktsqrindex << " / " << yindex << ". "
            << LINEINFO << endl;
        return;
    }
    /*if (n[ktsqrindex][yindex]>=-eps)
        cerr << "Overwriting data, ktindex = " << ktsqrindex << ", yindex = "
            << yindex << ". " << LINEINFO << endl; */
    n[ktsqrindex][yindex]=value;

    if (yindex>0)
        derivatives[ktsqrindex][yindex-1]=der;
    else
        cerr << "Added data point with yindex=" << yindex << " " << LINEINFO << endl;
}

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

    if (abs(ktsqr - par->ktsqr) < 1e-12)
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
            / abs(ktsqr - par->ktsqr);
    
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
    if (abs(ktsqr - par->ktsqr) < 1e-12)
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

REAL Amplitude::RapidityDerivative(REAL ktsqr, REAL y)
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

    REAL result, abserr; size_t eval;
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
void Amplitude::Solve(REAL maxy)
{
    // Find maxyind corresponding to maxy
    int maxyind=POINTS_Y-1;
    for (int i=1; i<POINTS_Y; i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }
    int largedifference=0;
    for (int yind=1; yind<maxyind; yind++)
    {
        cerr << "Solving for y=" << yvals[yind] << endl;
        // Solve N(y+DELTA_Y, kt) for every kt

#pragma omp parallel for
        for (int ktsqrind=0; ktsqrind<POINTS_KTSQR-1; ktsqrind++)
        {
            REAL tmpkt = ktsqrvals[ktsqrind];
            REAL tmpder = RapidityDerivative(tmpkt, yvals[yind-1]);
            REAL dy = yvals[yind]-yvals[yind-1];
            REAL newn=n[ktsqrind][yind-1] + dy*tmpder;

            
            
            AddDataPoint(ktsqrind, yind, newn, tmpder );
            if( abs(newn - n[ktsqrind][yind-1])/n[ktsqrind][yind-1] > 0.02)
            {
                largedifference++;
            }
            //cout << "N(ktsqr=" << ktsqrvals[ktsqrind] <<", y=" << yvals[yind] << ") = " << newn
            //<< ",  at lower rapidity it was " << n[ktsqrind][yind-1] << endl;
            
        }
    }
    cout << endl << "#" << largedifference << " out of " << maxyind * (POINTS_KTSQR-1)
            << " too large differences" << endl;
    // Again
    for (int avg=0; avg<1; avg++)
    {
        
        largedifference=0;
        for (int yind=1; yind<maxyind-1; yind++)
        {
            //cout << "Solving for y=" << yvals[yind] << endl;
            // Solve N(y+DELTA_Y, kt) for every kt
#pragma omp parallel for
            for (int ktsqrind=0; ktsqrind<POINTS_KTSQR-1; ktsqrind++)
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
        cout << endl << "#" << largedifference << " out of " << maxyind * (POINTS_KTSQR-1)
            << " too large differences" << endl;
    }

}


REAL Amplitude::InitialCondition(REAL ktsqr)
{
    switch(ic)
    {
        case INVPOWER:    
            return pow(ktsqr+1, -1);    // BK in Full Momentum Space, hep-ph/0504080
            break;
        case FTIPSAT:
            if (ktsqr > 2500) return 0.0;
            return gsl_sf_gamma_inc(0.0,ktsqr/4.0); // Ft of 2(1-exp(-r^2))
            break;
        default:
            cerr << "Unrecognized initial condition " << ic << endl;
    }
    //return exp(-0.5*SQR(log(ktsqr)-1)); // hep-ph/0110325 (bk+kinematical constraint)
}

void Amplitude::SetInitialCondition(INITIAL_CONDITION i)
{
    ic=i;
}

string Amplitude::InitialConditionStr()
{
    switch (ic)
    {
        case INVPOWER:
            return "(kt^2 + 1)^(-1), BK in full momentum space, hep-ph/0504080";
            break;
        case FTIPSAT:
            return "Gamma[0, kt^2/4], FT of 2(1-exp(-r^2))";
            break;
        default:
            return "Unknown initial condition";
            break;
    }
}

REAL Amplitude::Ktsqrval(int i)
{
    if (i > ktsqrvals.size()-1 or i<0)
    {
        cerr << "ktsqr index " << i <<" is out of range! " << LINEINFO << endl;
    }
    return ktsqrvals[i];
}

REAL Amplitude::Yval(int i)
{
    if (i > yvals.size()-1 or i<0)
    {
        cerr << "yval index " << i <<" is out of range! " << LINEINFO << endl;
    }
    return yvals[i];
}

void Amplitude::SetKinematicConstraint(bool kc)
{
    kinematic_constraint=kc;
}
