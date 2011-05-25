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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
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
            tmpvec.push_back(-1.0);
            tmpdervec.push_back(0.0);
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
    if (y<eps and datafile==false) return InitialCondition(ktsqr);
	if (y<0) y=0;
    // Linear interpolation or even extrapolation
    // TODO: SPLINE

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
        yind=yvals.size()-2;
        cerr << "Asked amplitude at too large Y=" << y << ", falling back to "
            << " y=" << yvals[yind] << ". " << LINEINFO << endl;
    }
    for (unsigned int i=0; i<=ktsqrvals.size()-2; i++)
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
		
		if (n[ktsqrind][yind+1]>-eps and y-yvals[yind]>eps)
        {	// Interpolate in y linearly
			tmparray[i-interpolation_start]=n[i][yind] 
			 + (y - yvals[yind]) * (n[i][yind+1] - n[i][yind]) / (yvals[yind+1]-yvals[yind]); 
		} 
		else
			tmparray[i-interpolation_start] = n[i][yind];
    }
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, interpo_points);
    gsl_spline_init(spline, tmpxarray, tmparray, interpo_points);

    REAL res;
    int status = gsl_spline_eval_e(spline, ktsqr, acc, &res);
    if (status)
    {
        cerr << "Interpolatioin failed at " << LINEINFO << ", error " << status << endl;
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    delete[] tmparray;
    delete[] tmpxarray;
    //return res;

    
    REAL linear = n[ktsqrind][yind] +
        (ktsqr - ktsqrvals[ktsqrind]) *
            (n[ktsqrind+1][yind] - n[ktsqrind][yind]) / (ktsqrvals[ktsqrind+1] - ktsqrvals[ktsqrind]);
    if (n[ktsqrind][yind+1]>-eps)   // We can also interpolate y
        linear += (y - yvals[yind]) * (n[ktsqrind][yind+1] - n[ktsqrind][yind]) / (yvals[yind+1]-yvals[yind]);
    

    if (std::abs(res-linear)/res > 0.5 and linear>eps)
		cout << "At ktsqr=" << ktsqr << ", y=" << y << " reldifference " << std::abs(res-linear)/res << endl;
    //return linear;
    
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
    /*if (n[ktsqrindex][yindex]>=-eps)
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
REAL Amplitude::LogLogDerivative(REAL ktsqr, REAL y)
{
    REAL epsilon = ktsqr/1000.0;
    return ktsqr/epsilon*log(N(ktsqr+epsilon,y)/N(ktsqr, y) );
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
	SetMaxKtsqr(minktsqr*pow(KtsqrMultiplier(), f.KtsqrPoints()));	
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
