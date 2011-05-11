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
    for (unsigned int i=0; i<KtsqrPoints(); i++)
        ktsqrvals.push_back(minktsqr * std::pow(ktsqr_multiplier, (int)i) );
    for (unsigned int i=0; i<YPoints()+1; i++)
        yvals.push_back((REAL) i * delta_y);
    
    for (unsigned int i=0; i<KtsqrPoints(); i++)   // Intialize every kt
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
		interpolation_end = KtsqrPoints()-2;
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
       // if (ktsqrind == 2734)
       //     cout << "Setting tmpxarray[" << i-interpolation_start << " to "
       // <<ktsqrvals[i] << ", i=" << i << " ktsqr=" << ktsqr << endl;
		tmpxarray[i-interpolation_start]=ktsqrvals[i];
		
		if (n[ktsqrind][yind+1]>-eps and y-yvals[yind]>eps)
        {	// Interpolate in y linearly
			tmparray[i-interpolation_start]=n[i][yind] 
			 + (y - yvals[yind]) * (n[i][yind+1] - n[i][yind]) / (yvals[yind+1]-yvals[yind]); 
		} 
		else
			tmparray[i-interpolation_start] = n[i][yind];
    }

    if (interpolation_start>2000 and tmpxarray[0]<1)
        cerr <<"tmpxarray[0]=" << tmpxarray[0] << " at y="<<y<< ", ktsqr=" << ktsqr
        << " ktsqrind=" << ktsqrind << ", interpolation_start: " << interpolation_start << endl;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, interpo_points);
    gsl_spline_init(spline, tmpxarray, tmparray, interpo_points);

    REAL res = gsl_spline_eval(spline, ktsqr, acc);
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
 * Interpolate
 * TODO, or remove completely?
 */
void Amplitude::Interpolate()
{
    

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
 * Solve BK using GSL ODE solver
 */
 

/*
 * Defines the equation for use in the ODE solver
 * t: "time"=rapidity, y[] = table of the amplitude values got in
 * previous iteration, f[] = table where to save new values of dy/dt = dN/dy,
 * params = other parameters
 */

struct inthelper_gslsolve
{
    REAL ktsqr;
    REAL* n;
};

REAL inthelperf_gslsolve( REAL ktsqr, void* p)
{

    return 0;
}

int func(double t, const double y[], double f[], void *params)
{
    double dummy1; 
    void * dummy2; 
    dummy1 =t;   // We don't need these
    dummy2 =params;

    inthelper_gslsolve inthelp;
    int *dim = (int*) params;
    inthelp.n = new REAL[(*dim)];
    for (int i=0; i< (*dim); i++)
        inthelp.n[i]=y[i];

    //Compute rapidity derivatives
    for (int i=0; i<(*dim); i++)
    {
        gsl_function fun;
        fun.params=&inthelp;
        fun.function = inthelperf_gslsolve;

        REAL result, abserr; 
        gsl_integration_workspace *workspace 
         = gsl_integration_workspace_alloc(KTSQRINTITERATIONS); 
        int status=gsl_integration_qag(&fun, 1e-8, 1e8, 0, KTSQRINTACCURACY, 
            KTSQRINTITERATIONS, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);

        f[i] = 0.2*(result - y[i]*y[i]);
    }

   /*
   for( i=0 ; i<=uvMax ; i++) {
      ws[i] = y[i];
   }
   for( i=0 ; i<=uvMax ; i++) {
      f[i] = abar[i] * cheb_approx(i);
   }*/
   return GSL_SUCCESS;
}


void Amplitude::SolveGSL(REAL maxy)
{
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, KtsqrPoints());
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.1, 0.0);
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (KtsqrPoints());

    int dim = KtsqrPoints();
    gsl_odeiv_system sys = {func, NULL, KtsqrPoints(), &dim};

    REAL* phi = new REAL[KtsqrPoints()];
    for (int i=0; i<KtsqrPoints(); i++)
    {
        phi[i]=n[i][0];
    }

    // Find maxyind corresponding to maxy
    int maxyind=YPoints();
    for (unsigned int i=1; i<=YPoints(); i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }

    for (int yind=1; yind <= maxyind; yind++)
    {
        REAL Y=0;
        REAL tmpy=yvals[yind];
        REAL h = yvals[yind]-yvals[yind-1];
        while( Y < tmpy)
        {
            int status = gsl_odeiv_evolve_apply(e, c, s, &sys,
                &Y, tmpy, &h, phi);
            if (status != GSL_SUCCESS)
            {
                cerr << "Error at " << LINEINFO << ": code " << status << endl;
            }
            yvals[yind]=Y;
            for (int i=0; i<KtsqrPoints(); i++)
                n[yind][i] = phi[i];

        }

    }

    delete[] phi;


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
