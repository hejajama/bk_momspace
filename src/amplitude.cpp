/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
//#include <bci.h>
#include <vector>
#include <sstream>
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
#include <gsl/gsl_roots.h>
#include <tools/interpolation.hpp>
#include "hankel.hpp"


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
    interpolation_rapidity=-1.0;
    interpolation_points = INTERPOLATION_POINTS;

    running_coupling = CONSTANT;
    alphas_scaling=1.0;
}

/*
 * Clear all
 */
void Amplitude::Clear()
{
	ktsqrvals.clear();
    ktvals.clear();
    lnktsqrvals.clear();
	yvals.clear();

    for (unsigned int i=0; i<n.size(); i++)
        n[i].clear();
    n.clear();
    derivatives.clear();
 
}

/*
 * Initialize
 * Initial condition must be set
 * 
 * NOTE: This function MUST BE called before any value is obtained from
 * this class.
 */
void Amplitude::Initialize()
{
	Clear();
    ktsqrpoints =  static_cast<unsigned int>(std::log(maxktsqr/minktsqr) 
            / std::log(ktsqr_multiplier) );
    maxktsqr = minktsqr * std::pow(ktsqr_multiplier, (int)ktsqrpoints);
    for (unsigned int i=0; i<KtsqrPoints(); i++)
    {
        ktsqrvals.push_back(minktsqr * std::pow(ktsqr_multiplier, (int)i) );
        lnktsqrvals.push_back(std::log(ktsqrvals[i]));
        ktvals.push_back(std::sqrt(ktsqrvals[i]));
    }

    yvals.push_back(0.0);

    // Intialize ln_n at y=0
    std::vector<REAL> tmpvec;
    std::vector<REAL> tmpdervec;
    for (uint i=0; i<KtsqrPoints(); i++)
    {
        REAL icval = InitialCondition( ktsqrvals[i] );
        tmpvec.push_back(icval);
        tmpdervec.push_back(0.0);
    }
    n.push_back(tmpvec);
    derivatives.push_back(tmpdervec);
}

/*
 * Compute amplitude from the tabulated values
 */

REAL Amplitude::N(REAL ktsqr, REAL y, bool bspline, bool derivative)
{
   
    if (ktsqr > 1.01*ktsqrvals[KtsqrPoints()-1] or 1.01*ktsqr < ktsqrvals[0])
        cerr << "Ktsqrval " << ktsqr << " is out of range [" << 
            MinKtsqr() << ", "  << MaxKtsqr() << "]! " << LINEINFO << endl;
    
    if (y<eps and datafile==false and derivative==false) return InitialCondition(ktsqr);
    if (y<eps) y=0;

    int ktsqrind;

    // Force ktsqr in range ]MinKtsqr(), MaxKtsqr()[
    
    if (ktsqr>=ktsqrvals[KtsqrPoints()-1]) ktsqr = ktsqrvals[KtsqrPoints()-1]*0.999;
    
    if (ktsqr <= ktsqrvals[0])
        ktsqr = ktsqrvals[0]*1.001;
        
    // Find ktsqrval and yval indexes refer to index for which
    // val[index]  is smaller than (or equal) y/ktsqr
    ktsqrind = FindIndex(ktsqr, ktsqrvals);

    REAL kt = std::sqrt(ktsqr);
    
    // Now we always have MinKtsqr() < ktsqr < MaxKtsqr() => interpolation works

    int yind = FindIndex(y, yvals);
    
    if (yind < 0 or ktsqrind<0) 
    {
        cerr << "Something very weird happened at " << LINEINFO << " with "
        << "ktsqr=" << ktsqr << ", y=" << y << endl;
        return 0;
    }
        
    // Keep y fixed, interpolate ktsqr
    // Interpolate only INTERPOLATION_POINTS points in order to make this
    // almost fast
    // Interpolate linearly in y and use spline in ktsqr

    int interpolation_start, interpolation_end;
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

    // First data point is sometimes somehow off, so don't use spline
    // with that
    ///TODO: Why?
    if (ktsqrind>1 and interpolation_start==0) { interpolation_start=1; }
    
	int interpo_points = interpolation_end - interpolation_start+1;
    
    REAL *tmparray = new REAL[interpo_points];
    REAL *tmpxarray = new REAL[interpo_points];
    for (int i=interpolation_start; i<= interpolation_end; i++)
    {
		tmpxarray[i-interpolation_start]=ktvals[i];

        REAL n0 = n[yind][i];
        tmparray[i-interpolation_start] = n0;

        // Interpolate in y if possible
        
		if (yind < yvals.size()-1 )
        {
            tmparray[i-interpolation_start]=n0
                + (y - yvals[yind]) * (n[yind+1][i] - n0)
                / (yvals[yind+1]-yvals[yind]);
		}
    }
    
    Interpolator interp(tmpxarray, tmparray, interpo_points);
    if (bspline)
        interp.SetMethod(INTERPOLATE_BSPLINE);
    interp.Initialize();
    REAL res;
    if (derivative)
    {
        res=interp.Derivative(kt)/(2.0*kt);
    }
    else
        res = interp.Evaluate(kt);

    delete[] tmparray;
    delete[] tmpxarray;
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

    if (ktsqrindex < 0 or yindex < 0 or ktsqrindex > n[yindex].size()-1)
    {
        cerr << "Invalid ktsqr/y index " << ktsqrindex << " / " << yindex << ". "
            << LINEINFO << endl;
        return;
    }
    n[yindex][ktsqrindex]=value;

    if (yindex>=1)
        derivatives[yindex-1][ktsqrindex]=der;
    else
        cerr << "Added data point with yindex=" << yindex << " " << LINEINFO << endl;
}


/*
 * d ln N(ktsqr) / d ln(ktsqr) = d ktsqr / d ln Ktsqr d ln N(ktsqr) /d ktsqr
 * = ktsqr d ln N(ktsqr) / d ktsqr = ktsqr/N * dN/dktsqr
 * =  ktsqr/eps ln [ N(ktsqr + eps) / N(ktsqr) ]
 */
struct Derivhelper_loglog
{
    REAL y;
    Amplitude* N;
};

REAL Derivhelperf_loglog(REAL ktsqr, void* p)
{
    Derivhelper_loglog* par = (Derivhelper_loglog*)p;
    //REAL ktsqr = std::exp(lnktsqr);
    return par->N->N(ktsqr, par->y, true) ;
}

REAL Amplitude::LogLogDerivative(REAL ktsqr, REAL rapidity)
{

    int ktsqrind = FindIndex(ktsqr, ktsqrvals);
    if (ktsqrind==0 or ktsqrind >= KtsqrPoints())
    {
        cerr << "Can't calculate derivative at the edge of the ktsqr range. "
                << LINEINFO << endl;
        return 0;
    }

    /*bool bspline=true;
    SetInterpolationPoints(20);
    return N(ktsqr, rapidity, bspline, true);
    return ktsqr/N(ktsqr, rapidity, bspline)*N(ktsqr, rapidity, bspline, true);
    */

    // f'(x) \approx [ f(x+h) - f(x-h) ] / 2h + o(h^2)
    //SetInterpolationPoints(INTERPOLATION_POINTS_DER);
    SetInterpolationPoints(150);

    // h=lnktsqr[ind+1]-lnktsqr[ind] = lnktsqr[ind] - lnktsqr[ind-1]
    REAL lnktsqr1 = std::log(ktsqrvals[ktsqrind-1]);
    REAL lnktsqr2 = std::log(ktsqrvals[ktsqrind+1]);
    REAL n1 = std::log(N(ktsqrvals[ktsqrind-1], rapidity, true) );
    REAL n2 = std::log(N(ktsqrvals[ktsqrind+1], rapidity, true) );
    REAL h = lnktsqr2 - lnktsqr1;
    REAL der = (n2-n1)/h;

    return der;
    //return ktsqr/N(ktsqr, rapidity)*der;

    /*
    
    const bool bspline = true;
    SetInterpolationPoints(INTERPOLATION_POINTS_DER);
    REAL n1 = N(ktsqr, rapidity, true);
    REAL n2;
    REAL ktsqr2;
    if (ktsqr + ktsqr/1000.0 < MaxKtsqr())
        ktsqr2 = ktsqr + ktsqr/1000.0;
    else
        ktsqr2 = ktsqr - ktsqr/1000.0;
	
	n2 = N(ktsqr2, rapidity, bspline);
    SetInterpolationPoints(INTERPOLATION_POINTS);

   return ktsqr/n1*(n2-n1)/(ktsqr2 - ktsqr);
    */
    /*
    // GSL
    //REAL lnktsqr = std::log(ktsqr);
    //REAL h = lnktsqr/1000.0;
    REAL h = (ktsqr - Ktsqrval( KtsqrIndex(ktsqr) + 1) )/2.0;

    REAL result,abserr;
    Derivhelper_loglog helper;
    helper.y=rapidity; helper.N=this;
    gsl_function fun;
    fun.function=&Derivhelperf_loglog;
    fun.params = &helper;
    gsl_deriv_central(&fun, ktsqr, h, &result, &abserr);

    if (std::abs(abserr/result)>0.2)
    {
        cerr << "Numerical derivation failed at " << LINEINFO
            << ": result=" << result << ", relerr=" << std::abs(abserr/result)
            << ", ktsqr=" << ktsqr << ", y=" << rapidity << endl;
    }
    
    return result*ktsqr/n1;;
    */

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

REAL SaturationHelperf(REAL lnktsqr, void* p)
{
   
    Sathelper* par = (Sathelper*)p;
    REAL val = -std::pow(std::exp(lnktsqr), par->gammac)
        *par->N->N(std::exp(lnktsqr), par->y /*, true */);
    return val;
    
}
REAL Amplitude::SaturationScale(REAL y)
{
    Sathelper helper;
    helper.y=y; helper.N=this;
    helper.gammac = 0.6275;
    helper.gammac = 0.5;
    
    gsl_function f; f.function=&SaturationHelperf;
    f.params=&helper;

    const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);

    // Interpolation doesn't work well at the edge of the ktsqr range, so limit
    // the studied interval
    REAL interval_min = Ktsqrval(1);
    REAL interval_max = Ktsqrval(KtsqrPoints()-2);
    REAL pos = SolveKtsqr(y, SATSCALE_N);

    // Find log of saturation scale
    pos = std::log(pos);
    interval_min = std::log(interval_min);
    interval_max = std::log(interval_max);

    gsl_error_handler_t* handler = gsl_set_error_handler_off();
    if (gsl_min_fminimizer_set(s, &f, pos, interval_min, interval_max)
        == GSL_EINVAL)
    {
        // Endpoints do not enclose a minimum (on the other hand
        // helper(pos) > helper(min),helper(max), but we can anyway continue
    }
    gsl_set_error_handler(handler);


    int status; int max_iter = 500; int iter=0;
    do
    {
        iter++;
        gsl_min_fminimizer_iterate(s);
        pos = gsl_min_fminimizer_x_minimum(s);
        interval_min = gsl_min_fminimizer_x_lower (s);
        interval_max = gsl_min_fminimizer_x_upper (s);
        ///TODO: relerror/abserror from ktsqrval difference
        status = gsl_min_test_interval (interval_min, interval_max, 0.0, 0.01);
        
    } while (status == GSL_CONTINUE and iter < max_iter);

    if (status == GSL_CONTINUE)
    {
        cerr << "Didn't find saturation scale when y=" << y  << ", iterated "
            << iter << "/" << max_iter << " times at " << LINEINFO << endl;
    }

    gsl_min_fminimizer_free(s);

    //cout << "start: " << std::exp(0.4*y) << " min: " << Ktsqrval(1) << " max: " << Ktsqrval(KtsqrPoints()-2)
    //    << " minimumsqr: " << pos << " iterations: " << iter << endl;

    return std::sqrt(std::exp(pos));
}

/*
 * Return ktsqr for which N(ktsqr, y)=amp
 * Used in another definition of the saturation scale
 */

struct KtsqrSolverHelper
{
    Amplitude *N;
    REAL y;
    REAL amp;
};

REAL KtsqrSolverHelperF(REAL x, void* p)
{
    KtsqrSolverHelper* par = (KtsqrSolverHelper*)p;
    return par->N->N(x, par->y)-par->amp;
}

// Find root N(ktsqr, y) - amp = 0
REAL Amplitude::SolveKtsqr(REAL y, REAL amp)
{
    const int MAX_ITER = 200;
    const REAL ROOTFINDACCURACY = 0.001;
    KtsqrSolverHelper help;
    help.N=this; help.y=y; help.amp=amp;
    gsl_function f;
    f.params = &help;
    f.function = &KtsqrSolverHelperF;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    
    gsl_root_fsolver_set(s, &f, ktsqrvals[0], ktsqrvals[ktsqrvals.size()-1]);
    /*cout << "y=" << y <<": min: " << KtsqrSolverHelperF(ktsqrvals[0], &help)
        << " max: " << KtsqrSolverHelperF(ktsqrvals[ktsqrvals.size()-1], &help)
        << " minktsqr: " << ktsqrvals[0] << " maxktsqr: " << ktsqrvals[ktsqrvals.size()-1] << endl;
    */
    int iter=0; int status; REAL min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ROOTFINDACCURACY);
        

    } while (status == GSL_CONTINUE and iter < MAX_ITER);

    if (iter>=MAX_ITER)
    {
        cerr << "Solving failed at y=" << y << " at " << LINEINFO << endl;
    }

    REAL res = gsl_root_fsolver_root(s);

    gsl_root_fsolver_free(s);

    return std::sqrt(res);

}

/*
 * Solve saturation scale in r space
 * Defined as N(r=1/Q_s) = const
 */
const REAL SATSCALE_R = 0.2;
struct SatscaleRHelper
{
    Hankel *transf;
    REAL y;
};

REAL SatscaleRHelperf(REAL lnr, void* p)
{
    SatscaleRHelper* par = (SatscaleRHelper*)p;
    return par->transf->Amplitude_r(std::exp(lnr), par->y) - SATSCALE_R;
}

REAL Amplitude::SaturationScaleR(REAL y)
{
    const int MAX_ITER = 30;
    const REAL ROOTFINDACCURACY = 0.05;
    
    Hankel transform(this);
    SatscaleRHelper help; help.transf = &transform;
    help.y=y;
    gsl_function f; f.params=&help;
    f.function=SatscaleRHelperf;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    REAL minlnr = std::log(1e-5);
    REAL maxlnr = std::log(10);
    gsl_root_fsolver_set(s, &f, minlnr, maxlnr);

    int iter=0; int status; REAL min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ROOTFINDACCURACY);
        //cout << "y=" << y <<", interval " << std::exp(min) << " - " << std::exp(max) << endl;

    } while (status == GSL_CONTINUE and iter < MAX_ITER);

    REAL res = gsl_root_fsolver_root(s);

    if (iter>=MAX_ITER)
    {
        cerr << "Solving failed at y=" << y << " at " << LINEINFO <<
            " result " << std::exp(res) << " relaccuracy "
            << (std::exp(max)-std::exp(min))/std::exp(res) << endl;
    }


    gsl_root_fsolver_free(s);

    return std::exp(-res);
    
}

REAL Amplitude::InitialCondition(REAL ktsqr)
{
    int status;
    x0=0.01;
    switch(ic)
    {
        case INVPOWER:    
            return pow(ktsqr+1, -1);    // BK in Full Momentum Space, hep-ph/0504080
            break;
        case FTIPSAT:
            // Ft of 1-exp(-r^2*Q_s0^2 / 4) = 1/2*Gamma[0, k^2/Q_s0^2]
            // Fitted to HERA data in arXiv:0902.1112: Q_s0^2 = 0.24GeV^2
            // Q_s0^2 defined in amplitude.hpp
            // NB: In ref. \alpha(s) depends on log(4C^2/r^2lambdaqcd^2), but
            // tools.cpp:Alpha_s depends on log(Q^2/lambdaqcd^2). Fitted value
            // for C^2 is 5.3
            
            gsl_sf_result res;
            status = gsl_sf_gamma_inc_e(0.0, ktsqr/Q0SQR, &res);
            if (status==15) return 0.0; // Overflow
            if (status)
                cerr << "gsl_sf_gamma_inc_e failed at " << LINEINFO << " with "
                    << "code " << status << " res " << 0.5*res.val << endl;
            return 0.5*res.val;
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
    std::stringstream s;
    switch (ic)
    {
        case INVPOWER:
            s << "(kt^2 + 1)^(-1), BK in full momentum space, hep-ph/0504080";
            break;
        case FTIPSAT:
            s << "0.5*Gamma[0, kt^2/Q_s0^2], FT of 1-exp(-r^2 Q_s0^2/4),"
                << " Q_s0^2 = " << Q0SQR << " GeV^2";
            break;
        case INVPOWER4:
            s << "(kt^4+1)^(-1), arbitrary";
            break;
        case GAUSS:
            s << "Exp[-(Log[k^2] + 2)^2/5], hep-ph/0110325";
            break;
        default:
            return "Unknown initial condition";
            break;
    }
    return s.str();
}

string Amplitude::RunningCouplingStr()
{
    if (datafile==true)
	{
		return "Data is read from a file, don't know what running coupling was used";
	}
    std::stringstream s;
    switch (running_coupling)
    {
    case CONSTANT:
        s << "Constant, ALPHABAR_s=" << ALPHABAR_s << " ";
        break;
    case PARENT_DIPOLE:
        s << "Parent dipole";
        break;
    case MINK:
        s << "MINK";
        break;
    case MAXK:
        s << "MAXK";
        break;
    }
    if (running_coupling!=CONSTANT)
        s << " \\alpha_s scaling factor: " << alphas_scaling;
    return s.str();
}


int Amplitude::ReadData(string file)
{
	Clear();
	DataFile f(file);
	SetMinKtsqr(f.MinKtsqr());
	SetKtsqrMultiplier(f.KtsqrMultiplier());
	SetMaxKtsqr(minktsqr*std::pow(KtsqrMultiplier(), f.KtsqrPoints()));
	Initialize();

	f.GetData(n, yvals);

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
    return yvals.size()-1;
}

unsigned int Amplitude::KtsqrPoints()
{
    return ktsqrpoints;
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
    return minktsqr*std::pow(ktsqr_multiplier, (int)(KtsqrPoints()-1));
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

void Amplitude::SetRunningCoupling(RUNNING_COUPLING rc)
{
    running_coupling=rc;
}

RUNNING_COUPLING Amplitude::RunningCoupling()
{
    return running_coupling;
}

void Amplitude::SetInterpolationPoints(int p)
{
    interpolation_points=p;
}

/*
 * Return (effective) rapidity from xbj
 * Depends on x_0 which is again determined by the initial condition
 */
REAL Amplitude::Y(REAL xbj)
{
    return std::log(x0/xbj);
}

/* Returns index i for which
 * vec[i]<=val
 * If such index can't be found, returns -1
 */

int Amplitude::FindIndex(REAL val, std::vector<REAL> &vec)
{
    if (val < vec[0]) return -1;
    
    int ind=-1;
    
    uint start=0; uint end=vec.size()-1;
    while(end-start>10)
    {
        int tmp = static_cast<int>((start+end)/2.0);
        
        if (vec[tmp]>=val)
            end=tmp;
        else
            start=tmp;
    }
    
    
    for (uint i=start; i<=end; i++)
    {
        if (vec[i]<=val and vec[i+1]>val)
        {
            ind=i;
            break;
        }
    }
    if (ind == -1) return vec.size()-1;
    return ind;
}




/*
 * Add new rapidity value for tables
 */
void Amplitude::AddRapidity(REAL y)
{
    if (y <= yvals[yvals.size()-1])
    {
        cerr << "Trying to add rapidity value " << y << " but currently "
        << "the largest rapidity tabulated is " << yvals[yvals.size()-1] << endl;
        return;
    }

    yvals.push_back(y);

    std::vector<REAL> tmpvec;
    std::vector<REAL> tmpdervec;
    for (uint i=0; i<KtsqrPoints(); i++)
    {
        tmpvec.push_back(0.0);
        tmpdervec.push_back(0.0);
    }
    n.push_back(tmpvec);
    derivatives.push_back(tmpdervec);

}

void Amplitude::SetAlphasScaling(REAL scaling)
{
    alphas_scaling = scaling;
}

REAL Amplitude::AlphasScaling()
{
    return alphas_scaling;
}
