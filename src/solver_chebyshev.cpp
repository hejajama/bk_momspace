/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "solver_chebyshev.hpp"
#include <cmath>
#include <gsl/gsl_math.h>
using std::cos;
using std::log;
using std::exp;

void ChebyshevSolver::Solve(REAL maxy)
{
    Prepare();
    
    // Find maxyind corresponding to maxy
    unsigned int maxyind=YPoints();
    for (unsigned int i=1; i<=YPoints(); i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }

    for (unsigned int yind=1; yind <= maxyind; yind++)
    {
        REAL tmpy = yvals[yind];
        cout << "Solving for y=" << tmpy << endl;
        REAL dy = yvals[yind]-yvals[yind-1];
#pragma omp parallel for
        for (unsigned int uind=0; uind<KtsqrPoints()-1; uind++)
        {
            REAL derivative=0;
            derivative = Integrate(uind, yind-1);
            derivative -= SQR(N(ktsqrvals[uind], yvals[yind-1]));
            derivative *= 0.2;  ///TODO

            if (derivative < 0)
                cerr << "der: " << derivative << " at y=" << yind << ", uind=" << uind
                    <<", ktsqr=" << ktsqrvals[uind] << endl;
           

            REAL newn = N(ktsqrvals[uind], yvals[yind-1]) + dy*derivative;

            AddDataPoint(uind, yind, newn, derivative);

        }

    }


}


/*
 * Prepare all tables etc. for solution process
 */
void ChebyshevSolver::Prepare()
{
    uvals.clear();
    maxl = log(MaxKtsqr());
    minl = log(MinKtsqr());
    
    for (unsigned int i=0; i<=KtsqrPoints(); i++)
    {
        // NB: Small difference when compared with BKSolver
        uvals.push_back(cos(M_PI*(static_cast<REAL>(KtsqrPoints()-i)+0.5)/(2.0*KtsqrPoints())));
        if (i<KtsqrPoints())
            ktsqrvals[i]=exp( (-minl + maxl)*uvals[i]+minl );
        else
            ktsqrvals.push_back( (-minl + maxl)*uvals[i]+minl );
    }
}

/*
 * Compute integral appearing in \partial Y = \alpha_s \int [] - \alpha_s N^2
 * y = rapidity at which the derivative is computed
 */
REAL ChebyshevSolver::Integrate(int uind, int yind)
{
    // Values of the integrand
    REAL * f = new REAL[KtsqrPoints()+1];
    REAL result=0;
    
    for (unsigned int nind=0; nind<=CHEBYSHEV_DEGREE; nind++)
    {
        // Compute coefficients a_k
        ///TODO: USE FFT; KtsqrPoints() should be power of 2
        REAL a=0;
        for (unsigned int tmpk=0; tmpk < KtsqrPoints()-1; tmpk++)
        {
            a+=cos(nind*M_PI*(tmpk + 0.5)/(2.0*KtsqrPoints()))
                * Integrand( uind, tmpk, yind );
        }
        a *= 2.0/(2.0*KtsqrPoints());
        if (nind==0) a/=2.0;
        result += a * ChebyshevIntegral(nind);
    }  
    
    delete[] f;

    return result;
}

///TODO: Kinematical constraint
REAL ChebyshevSolver::Integrand(int uind, int vind, int yind)
{
    REAL result = 0;
    REAL l1 =  (-minl + maxl)*uvals[uind]+minl;  // Fixed
    REAL l2 =  (-minl + maxl)*uvals[vind]+minl;  // Integrated over
    if (uind != vind)
    {
        result += 1.0 / std::abs( 1.0 - exp(l1-l2) )
            * ( N(ktsqrvals[vind], yvals[yind]) - exp(l1-l2)*N(ktsqrvals[uind], yvals[yind]) );
    }
    result += 1.0 / sqrt(1.0 + 4.0*exp(l1-l2) ) * N(ktsqrvals[uind], yvals[yind]);
        
    return result;

}

/*
 * Compute \int_0^1 dx T_n(x)
 */
REAL ChebyshevSolver::ChebyshevIntegral(int n)
{
    if (n==1)
        return 0.5;
    else if (GSL_IS_ODD(n) == 1)
        return (-1.0 + n*sin(n*M_PI_2))/(-1.0+n*n);
    else
        return 1.0/(1.0-n*n);

}
