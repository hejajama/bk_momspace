/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
//#include <bci.h>
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_min.h>


#include <cmath>


using std::abs;

Amplitude::Amplitude()
{
    maxktsqr = DEFAULT_MAXKTSQR;
    minktsqr = DEFAULT_MINKTSQR;
    ktsqr_multiplier = DEFAULT_KTSQR_MULTIPLIER;
    maxy = DEFAULT_MAXY;
    delta_y = DEFAULT_DELTA_Y;
    averages=0;
    datafile=false;
    adams_method=false;
    interpolation_end=-1;
    interpolation_start=-1;
    interpolation_rapidity=-1.0;

    running_coupling = false;
}

/*
 * Clear all
 */
void Amplitude::Clear()
{
	ktsqrvals.clear();
	yvals.clear();

    for (unsigned int i=0; i<n.size(); i++)
        n[i].clear();
    n.clear();
 
}

/*
 * Initialize
 * Initial condition must be set
 */
void Amplitude::Initialize()
{
	Clear();
    for (unsigned int i=0; i<=KtsqrPoints(); i++)
        ktsqrvals.push_back(minktsqr * std::pow(ktsqr_multiplier, (int)i) );
    for (unsigned int i=0; i<YPoints()+1; i++)
        yvals.push_back((REAL) i * delta_y);
    
    for (unsigned int i=0; i<=KtsqrPoints(); i++)   // Intialize every kt
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
}

/*
 * Compute amplitude from the tabulated values
 */

REAL Amplitude::N(REAL ktsqr, REAL y, bool bspline)
{
    
    if (y<eps and datafile==false) return InitialCondition(ktsqr);
    if (y<eps) y=0;
    
    // Find ktsqrval and yval indexes refer to index for which
    // val[index]  is smaller than (or equal) y/ktsqr
    int yind = -1;
    int ktsqrind=-1;
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
    for (unsigned int i=0; i<ktsqrvals.size()-1; i++)
    {
        if (ktsqrvals[i]<=ktsqr and ktsqrvals[i+1]>ktsqr)
        {
            ktsqrind=i;
            break;
        }
    }
    if (ktsqrind < 0) // Didn't find, so refers to the largest one 
    {
        ktsqrind=ktsqrvals.size()-1;
        if (ktsqr - ktsqrvals[ktsqrvals.size()-1] > 100)
            cerr << "Asked amplitude at too large ktsqr=" << ktsqr << ", falling back to "
                << " ktsqr=" << ktsqrvals[ktsqrind] << ". " << LINEINFO << endl;
    }

    if (std::abs(ktsqr - ktsqrvals[0])/ktsqrvals[0] < 0.001)    // Don't interpolate in kt
    {
        if (yind == yvals.size()-1) return n[0][yind];
        return n[0][yind]+(y-yvals[yind])*(n[0][yind+1]-n[0][yind])
            / ( yvals[yind+1] - yvals[yind] );
    }

    // Keep y fixed, interpolate ktsqr
    // Interpolate only INTERPOLATION_POINTS points in order to make this
    // almost fast
    // Interpolate linearly in y and use spline in ktsqr

    unsigned int interpolation_start, interpolation_end;
    if (ktsqrind - INTERPOLATION_POINTS/2 < 0)
    {
		interpolation_start=0;
		interpolation_end=INTERPOLATION_POINTS;
	}
	else if (ktsqrind + INTERPOLATION_POINTS/2 > KtsqrPoints()-2 )
	{
		interpolation_end = KtsqrPoints()-1;
		interpolation_start = KtsqrPoints()-INTERPOLATION_POINTS-2;
	}
	else
	{
		interpolation_start = ktsqrind - INTERPOLATION_POINTS/2;
		interpolation_end = ktsqrind + INTERPOLATION_POINTS/2;
	}
	int interpo_points = interpolation_end - interpolation_start+1;
    
    REAL *tmparray = new REAL[interpo_points];
    REAL *tmpxarray = new REAL[interpo_points];
    for (int i=interpolation_start; i<= interpolation_end; i++)
    {
		tmpxarray[i-interpolation_start]=ktsqrvals[i];

        tmparray[i-interpolation_start] = n[i][yind];

        // Interpolate in y if possible
		if (yind < yvals.size()-1 )
        {
            if (n[ktsqrind][yind+1]>0 and y-yvals[yind]>0.000001 )
            {
                tmparray[i-interpolation_start]=n[i][yind] 
                 + (y - yvals[yind]) * (n[i][yind+1] - n[i][yind]) / (yvals[yind+1]-yvals[yind]);
            }
		} 
		
			
    }

    /*if (bspline==true)
    {
        REAL n = BSplineAmplitude(ktsqr, tmparray, tmpxarray, interpo_points);
        return n;
    }*/
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, interpo_points);
    gsl_spline_init(spline, tmpxarray, tmparray, interpo_points);

    REAL res;
    int status = gsl_spline_eval_e(spline, ktsqr, acc, &res);
    if (status)
    {
        cerr << "Interpolatioin failed at " << LINEINFO << ", error " << gsl_strerror(status)
         << " (" << status << "), ktsqr=" << ktsqr << ", y=" << y << endl;
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    delete[] tmparray;
    delete[] tmpxarray;
    //return res;
    /*
    REAL linear = n[ktsqrind][yind] +
        (ktsqr - ktsqrvals[ktsqrind]) *
            (n[ktsqrind+1][yind] - n[ktsqrind][yind]) / (ktsqrvals[ktsqrind+1] - ktsqrvals[ktsqrind]);
    if (n[ktsqrind][yind+1]>-eps)   // We can also interpolate y
        linear += (y - yvals[yind]) * (n[ktsqrind][yind+1] - n[ktsqrind][yind]) / (yvals[yind+1]-yvals[yind]);
    

    if (std::abs(res-linear)/res > 0.5 and linear>eps)
		cerr << "At ktsqr=" << ktsqr << ", y=" << y << " reldifference " << std::abs(res-linear)/res << endl;
    return linear;
    */
    return res;
    

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
    /*if (n[ktsqrindex][yindex]>    REAL n1,n2;
    
        cerr << "Overwriting data, ktindex = " << ktsqrindex << ", yindex = "
            << yindex << ". " << LINEINFO << endl; */
    n[ktsqrindex][yindex]=value;

    if (yindex>=0)
        derivatives[ktsqrindex][yindex-1]=der;
    else
        cerr << "Added data point with yindex=" << yindex << " " << LINEINFO << endl;
}


/*
 * d ln N(ktsqr) / d ln(ktsqr) = d ktsqr / d ln Ktsqr d ln N(ktsqr) /d ktsqr
 * = ktsqr d ln N(ktsqr) / d ksqr = ktsqr/eps ln [ N(ktsqr + eps) / N(ktsqr) ]
 */
struct Derivhelper_loglog
{
    REAL y;
    Amplitude* N;
};

REAL Derivhelperf_loglog(REAL lnktsqr, void* p)
{
    Derivhelper_loglog* par = (Derivhelper_loglog*)p;
    REAL ktsqr = std::exp(lnktsqr);
    return std::log(par->N->N(ktsqr, par->y) );
}

REAL Amplitude::LogLogDerivative(REAL ktsqr, REAL rapidity)
{

    //if (std::abs(rapidity-bspline_y)>0.001)
    //    IntializeBSpline(rapidity);

    REAL n1,n2;
    n1 = BSplineAmplitude(ktsqr, rapidity);
    n2 = BSplineAmplitude(ktsqr+ktsqr/100.0, rapidity);
	//n1 = N(ktsqr, rapidity);
	//n2 = N(ktsqr+ktsqr/100.0, rapidity);
    

   return ktsqr/n1*(n2-n1)/(ktsqr/100.0);
    
    /* Lowest order, not very accurate 
    REAL epsilon = ktsqr/1000.0;
    return ktsqr/epsilon*log(N(ktsqr+epsilon,y)/N(ktsqr, y) );
    */


    /* GSL
    REAL lnktsqr = std::log(ktsqr);
    REAL h = lnktsqr/100.0;

    REAL result,abserr;
    Derivhelper_loglog helper;
    helper.y=y; helper.N=this;
    gsl_function fun;
    fun.function=&Derivhelperf_loglog;
    fun.params = &helper;
    gsl_deriv_central(&fun, lnktsqr, h, &result, &abserr);

    if (std::abs(abserr/result)>0.1)
    {
        cerr << "Numerical derivation failed at " << LINEINFO
            << ": result=" << result << ", relerr=" << std::abs(abserr/result)
            << ", ktsqr=" << ktsqr << ", y=" << y << endl;
    }
    
    return result;
    */   

}

void Amplitude::IntializeBSpline(int ktsqrind, REAL rapidity)
{
    // Not first time -> free old memory
    if (interpolation_rapidity>=0.0)
    {
        gsl_bspline_free(bw);
        gsl_vector_free(B);
        gsl_matrix_free(X);
        gsl_vector_free(c);
        
        gsl_matrix_free(cov);
        gsl_multifit_linear_free(mw);
    }
    
    // Find ktsqrval and yval indexes refer to index for which
    // val[index]  is smaller than (or equal) y/ktsqr
    int yind = -1;
    for (unsigned int i=0; i<yvals.size()-1; i++)
    {
        if (yvals[i]<=rapidity and yvals[i+1]>rapidity)
        {
            yind=i;
            break;
        }
    }
    if (yind < 0) // Didn't find, so refers to the largest one 
    {
        yind=yvals.size()-1;
        if (rapidity-yvals[yind]>0.05)
            cerr << "Asked derivative at too large Y=" << rapidity << ", falling back to "
                << " y=" << yvals[yind] << ". " << LINEINFO << endl;
    }
    
    if (ktsqrind - INTERPOLATION_POINTS_DER/2 < 0)
    {
		interpolation_start=0;
		interpolation_end=INTERPOLATION_POINTS_DER;
	}
	else if (ktsqrind + INTERPOLATION_POINTS_DER/2 > KtsqrPoints()-2 )
	{
		interpolation_end = KtsqrPoints()-1;
		interpolation_start = KtsqrPoints()-INTERPOLATION_POINTS_DER-2;
	}
	else
	{
		interpolation_start = ktsqrind - INTERPOLATION_POINTS_DER/2;
		interpolation_end = ktsqrind + INTERPOLATION_POINTS_DER/2;
	}
	int interpo_points = interpolation_end - interpolation_start+1;


    gsl_vector *x, *y, *w;
    x = gsl_vector_alloc(interpo_points);
    y = gsl_vector_alloc(interpo_points);
    w = gsl_vector_alloc(interpo_points);

    for (int i=interpolation_start; i<= interpolation_end; i++)
    {
		//tmpxarray[i-interpolation_start]=ktsqrvals[i];
        gsl_vector_set(x, i-interpolation_start, ktsqrvals[i]);
        
        gsl_vector_set(w, i-interpolation_start, 1.0);

        //tmpxarray[i-interpolation_start]=ktsqrvals[i];
        REAL tmpy;
		
		if (n[ktsqrind][yind+1]>-eps and rapidity-yvals[yind]>eps)
        {	// Interpolate in y linearly
			tmpy=n[i][yind] 
			 + (rapidity - yvals[yind]) * (n[i][yind+1] - n[i][yind]) / (yvals[yind+1]-yvals[yind]); 
		}
        else
            tmpy = n[i][yind];
        gsl_vector_set(y, i-interpolation_start, tmpy);
    }
    
    const int ncoeffs = 12;
    const int nbreak = ncoeffs-2; // k=4
     
    /* allocate a cubic bspline workspace (k = 4) */
    bw = gsl_bspline_alloc(4, nbreak);
    B = gsl_vector_alloc(ncoeffs);
       
    X = gsl_matrix_alloc(interpo_points, ncoeffs);
    c = gsl_vector_alloc(ncoeffs);
       
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    mw = gsl_multifit_linear_alloc(interpo_points, ncoeffs);
     
     
    // use uniform breakpoints on 
    gsl_bspline_knots_uniform(ktsqrvals[interpolation_start],
            ktsqrvals[interpolation_end], bw);
     
    /* construct the fit matrix X */
    for (int i = 0; i < interpo_points; ++i)
    {
       double xi = gsl_vector_get(x, i);
     
       /* compute B_j(xi) for all j */
       gsl_bspline_eval(xi, B, bw);
     
       /* fill in row i of X */
       for (int j = 0; j < ncoeffs; ++j)
       {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(X, i, j, Bj);
       }
    }
     
    /* do the fit */
    REAL chisq;
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_vector_free(w);

    interpolation_rapidity=rapidity;
}



REAL Amplitude::BSplineAmplitude(REAL ktsqr, REAL rapidity)
{
    /*REAL n,yerr;
    gsl_bspline_eval(ktsqr, B, bw);
    gsl_multifit_linear_est(B, c, cov, &n, &yerr);
    return n;
    */

    
    // Find ktsqrval and yval indexes refer to index for which
    // val[index]  is smaller than (or equal) y/ktsqr
    int yind = -1;
    for (unsigned int i=0; i<yvals.size()-1; i++)
    {
        if (yvals[i]<=rapidity and yvals[i+1]>rapidity)
        {
            yind=i;
            break;
        }
    }
    if (yind < 0) // Didn't find, so refers to the largest one 
    {
        yind=yvals.size()-1;
        if (rapidity-yvals[yind]>0.05)
            cerr << "Asked derivative at too large Y=" << rapidity << ", falling back to "
                << " y=" << yvals[yind] << ". " << LINEINFO << endl;
    }
    int ktsqrind=-1;
    for (unsigned int i=0; i<ktsqrvals.size()-1; i++)
    {
        if (ktsqrvals[i]<=ktsqr and ktsqrvals[i+1]>ktsqr)
        {
            ktsqrind=i;
            break;
        }
    }
    if (ktsqrind < 0) // Didn't find, so refers to the largest one 
    {
        ktsqrind=ktsqrvals.size()-1;
        if (ktsqr - ktsqrvals[ktsqrvals.size()-1] > 100)
            cerr << "Asked derivative at too large ktsqr=" << ktsqr << ", falling back to "
                << " ktsqr=" << ktsqrvals[ktsqrind] << ". " << LINEINFO << endl;
    }

    /*if (ktsqrind < interpolation_start+10 or ktsqrind > interpolation_end - 10
        or std::abs(rapidity < interpolation_rapidity) > 0.01)
    {
        IntializeBSpline(ktsqrind, rapidity);
    }*/
    IntializeBSpline(ktsqrind, rapidity);
   
    

    REAL n1,yerr;
    gsl_bspline_eval(ktsqr, B, bw);
    gsl_multifit_linear_est(B, c, cov, &n1, &yerr);

    return n1;

}

/*
 * Saturation Scale
 * Calculates Q_s
 * Note: There area different definitions for Q_s, here the definition similar than
 * in hep-ph/0504080 is used:
 * k^(2*\gamma_c)N(k) has a maximum at k=Q_s, \gamma_c is the anomalous
 * dimension \gamma_c=0.6275
 * Note: in hep-ph/0504080 a reduced wave front with \gamma_c=0.5 is studied
 * due to the small lattice size, which is not a problem here if Y is small enough
 * (Y=40 is not too large)
 */

struct Sathelper
{
    Amplitude* N;
    REAL y;
    REAL gammac;
};

REAL SaturationHelperf(REAL ktsqr, void* p)
{
    Sathelper* par = (Sathelper*)p;
    //return -std::pow(ktsqr, par->gammac)*par->N->N(ktsqr, par->y);
    return -std::pow(ktsqr, par->gammac)*par->N->BSplineAmplitude(ktsqr, par->y);
    
}
REAL Amplitude::SaturationScale(REAL y)
{
    Sathelper helper;
    helper.y=y; helper.N=this;
    helper.gammac = 0.6275;
    
    gsl_function f; f.function=&SaturationHelperf;
    f.params=&helper;

    const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);

    // Interpolation doesn't work well at the edge of the ktsqr range, so limit
    // the studied interval
    REAL interval_min = Ktsqrval(1);
    REAL interval_max = Ktsqrval(KtsqrPoints()-2);
    REAL pos = std::exp( 0.4*y );   // ref: Q_s^2 \sim exp(v_c*Y)
    if (pos < interval_min) pos = 2.0*interval_min;
    if (pos>interval_max) pos = 0.5*interval_max;
    
    gsl_min_fminimizer_set(s, &f, pos, interval_min, interval_max);
    int status; int max_iter = 500; int iter=0;
    do
    {
        iter++;
        gsl_min_fminimizer_iterate(s);
        pos = gsl_min_fminimizer_x_minimum(s);
        interval_min = gsl_min_fminimizer_x_lower (s);
        interval_max = gsl_min_fminimizer_x_upper (s);
        status = gsl_min_test_interval (interval_min, interval_max, 0.0, 0.001);
        
    } while (status == GSL_CONTINUE and iter < max_iter);

    if (status == GSL_CONTINUE)
    {
        cerr << "Didn't find saturation scale when y=" << y  << ", iterated "
            << iter << "/" << max_iter << " times at " << LINEINFO << endl;
    }

    gsl_min_fminimizer_free(s);

    //cout << "start: " << std::exp(0.4*y) << " min: " << Ktsqrval(1) << " max: " << Ktsqrval(KtsqrPoints()-2)
    //    << " minimumsqr: " << pos << " iterations: " << iter << endl;

    return std::sqrt(pos);
}

/*
REAL Amplitude::BSplineDerivative(REAL ktsqr, REAL* ktsqrarray, REAL* narray, uint points)
{
    REAL n1,n2;
    n1 = BSplineAmplitude(ktsqr, ktsqrarray, narray, points);
    n2 = BSplineAmplitude(ktsqr+ktsqr/100.0, ktsqrarray, narray, points);
    

   return ktsqr/n1*(n2-n1)/(ktsqr/100.0);

}*/


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
        case INVPOWER4:
            return pow(SQR(ktsqr)+1, -1);
            break;
        case GAUSS:
            return std::exp( -SQR((std::log(ktsqr)+2))/5 );
            break;
        default:
            cerr << "Unrecognized initial condition " << ic << endl;
    }
    //return exp(-0.5*SQR(log(ktsqr)-1)); // hep-ph/0110325 (bk+kinematical constraint)
    return -1.0;
}

void Amplitude::SetInitialCondition(INITIAL_CONDITION i)
{
    ic=i;
}

string Amplitude::InitialConditionStr()
{
	if (datafile==true)
	{
		return "Data is read from a file, don't know what initial condition was used";
	}
    switch (ic)
    {
        case INVPOWER:
            return "(kt^2 + 1)^(-1), BK in full momentum space, hep-ph/0504080";
            break;
        case FTIPSAT:
            return "Gamma[0, kt^2/4], FT of 2(1-exp(-r^2))";
            break;
        case INVPOWER4:
            return "(kt^4+1)^(-1), arbitrary";
            break;
        case GAUSS:
            return "Exp[-(Log[k^2] + 2)^2/5], hep-ph/0110325";
        default:
            return "Unknown initial condition";
            break;
    }
}


int Amplitude::ReadData(string file)
{
	Clear();
	DataFile f(file);
	SetMinKtsqr(f.MinKtsqr());
	SetKtsqrMultiplier(f.KtsqrMultiplier());
	SetMaxKtsqr(minktsqr*std::pow(KtsqrMultiplier(), f.KtsqrPoints()));
	SetMaxY(f.MaxY());
	SetDeltaY(f.DeltaY());
	Initialize();

	f.GetData(n);

	datafile=true;
	return 0;
}


REAL Amplitude::Ktsqrval(unsigned int i)
{
    if (i > ktsqrvals.size()-1 or i<0)
    {
        cerr << "ktsqr index " << i <<" is out of range! " << LINEINFO << endl;
    }
    return ktsqrvals[i];
}

REAL Amplitude::Yval(unsigned int i)
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

void Amplitude::SetNumberOfAveragements(int avg)
{
    averages=avg;
}

void Amplitude::SetMinKtsqr(REAL mkt)
{
    minktsqr = mkt;
}

void Amplitude::SetMaxKtsqr(REAL mkt)
{
    maxktsqr = mkt;
}

void Amplitude::SetKtsqrMultiplier(REAL m)
{
    ktsqr_multiplier=m;
}

void Amplitude::SetMaxY(REAL y)
{
    maxy=y;

}

void Amplitude::SetDeltaY(REAL dy)
{
    delta_y = dy;
}

unsigned int Amplitude::YPoints()
{
    return static_cast<unsigned int>(maxy/delta_y)+1;
}

unsigned int Amplitude::KtsqrPoints()
{
    return static_cast<unsigned int>(std::log(maxktsqr/minktsqr) / std::log(ktsqr_multiplier) );
}

REAL Amplitude::DeltaY()
{
    return delta_y;
}

REAL Amplitude::KtsqrMultiplier()
{
    return ktsqr_multiplier;
}

REAL Amplitude::MinKtsqr()
{
    return minktsqr;
}

REAL Amplitude::MaxKtsqr()
{
    return maxktsqr;
}

bool Amplitude::KinematicalConstraint()
{
    return kinematic_constraint;
}

Amplitude::~Amplitude()
{

    // If bspline interpolation was used
    if (interpolation_rapidity>=0)
    {
        gsl_bspline_free(bw);
        gsl_vector_free(B);
        gsl_matrix_free(X);
        gsl_vector_free(c);
        
        gsl_matrix_free(cov);
        gsl_multifit_linear_free(mw);
    }
    
}

void Amplitude::SetRunningCoupling(bool rc)
{
    running_coupling=rc;
}

bool Amplitude::RunningCoupling()
{
    return running_coupling;
}
