/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "solver_chebyshev.hpp"
#include <cmath>
#include <gsl/gsl_math.h>


// Fourier coisne transofrmation from src/fft4g.c
extern "C"
{
    void   dfct  (int n, double *a, double *t, int *ip, double *w);
}

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
        for (unsigned int uind=0; uind<=KtsqrPoints(); uind++)
        {
            REAL derivative=0;
            derivative = Integrate(uind, yind-1);
            derivative -= SQR(N(ktsqrvals[uind], yvals[yind-1]));
            derivative *= 0.2;  ///TODO

            if (derivative < 0)
                cerr << "der: " << derivative << " at y=" << yind << ", uind=" << uind
                    <<", ktsqr=" << ktsqrvals[uind] << endl;

            REAL newn=0;
            // Adams method: apprximate derivative as a 2nd order polynomial
            // Y_{i+1} = y_i + hf_i + 1/2 h (f_i - f_{i-1} )
            // We can use this only for yind>1
            if (yind>1 and adams_method==true)
            {
                REAL old_der = derivatives[uind][yind-2];
                newn = N(ktsqrvals[uind], yvals[yind-1]) + dy*derivative
                    + 1.0/2.0*dy*( derivative - old_der);
            }
           else
                newn = N(ktsqrvals[uind], yvals[yind-1]) + dy*derivative;

            AddDataPoint(uind, yind, newn, derivative);

        }

    }


}


/*
 * Compute integral appearing in \partial Y = \alpha_s \int [] - \alpha_s N^2
 * y = rapidity at which the derivative is computed
 */
REAL ChebyshevSolver::Integrate(int uind, int yind)
{
    // Values of the integrand
    REAL * f = new REAL[2*KtsqrPoints()+1];
    

    // Compute coefficients a_k

    /* Fill half the array with the integrand,
      since we only use half the interval - could probably be done better
       KtsqrPoints() must be a power of 2
    */
    for(unsigned k=1 ; k<=KtsqrPoints() ; k++ ) {
        f[k] = Integrand(uind, KtsqrPoints()-k, yind);
        f[2*KtsqrPoints()-k] = 0.0;
    }
    f[KtsqrPoints()]=0;
    f[2*KtsqrPoints()]=0;

    //cout << f[0] << " " << f[1] << " " << f[KtsqrPoints()] << " " <<
    //f[KtsqrPoints()*2-2] << endl;
       
    // Call the FFT cosine transform
    dfct(2*KtsqrPoints(), f, twork, ipwork, wwork);

    REAL sum = 0.5*f[0]*ChebyshevIntegral(0);
    for (unsigned int k=1; k<2*KtsqrPoints(); k++)
        sum += f[k]*ChebyshevIntegral(k);

    sum *= 2.0/static_cast<REAL>(2*KtsqrPoints());

    sum *= (-minl + maxl);

    delete[] f;
    return sum;
        /*
        
        for (unsigned int nind=0; nind<=CHEBYSHEV_DEGREE; nind++)
        {
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
        */
    //}  
    
    //delete[] f;

   // return result;
}

///TODO: Kinematical constraint
REAL ChebyshevSolver::Integrand(int uind, int vind, int yind)
{
    REAL result = 0;
    REAL l1 =  (-minl + maxl)*uvals[uind]+minl;  // Fixed
    REAL l2 =  (-minl + maxl)*uvals[vind]+minl;  // Integrated over
    if (uind != vind)
    {
        result += 1.0 / std::abs( 1.0 - std::exp(l1-l2) )
            * ( N(ktsqrvals[vind], yvals[yind]) - std::exp(l1-l2)*N(ktsqrvals[uind], yvals[yind]) );
    }
    result += 1.0 / sqrt(1.0 + 4.0*std::exp(2.0*(l2-l1) ) ) * N(ktsqrvals[uind], yvals[yind]);
        
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
        return (-1.0 + n*std::sin(n*M_PI_2))/(-1.0+n*n);
    else
        return 1.0/(1.0-n*n);

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
        //uvals.push_back(cos(M_PI*(static_cast<REAL>(KtsqrPoints()-i)+0.5)/(2.0*KtsqrPoints())));
        uvals.push_back(std::cos(M_PI*(static_cast<REAL>(KtsqrPoints()-i))/(2.0*KtsqrPoints())));
      //  if (i < KtsqrPoints())
        ktsqrvals[i]=(std::exp( (-minl + maxl)*uvals[i]+minl ) );
        n[i][0] = InitialCondition(ktsqrvals[i]);
    //    else
    //        ktsqrvals.push_back( exp( (-minl + maxl)*uvals[i]+minl ) );
    }
}

/*
 * Constructor and destructor
 * Allocate memory
 */
ChebyshevSolver::ChebyshevSolver()
{
    wwork = new REAL[2*KtsqrPoints()+1];
    ipwork  = new int[2+(int)sqrt((2*KtsqrPoints()+1)/4.0)];
    twork = new REAL[2*KtsqrPoints()+1];

    // There must be 2^n ktsqrvals -> force this
    ///TODO: THIS IS STUPID
    minktsqr=1e-6;
    maxktsqr=1e6;
    int N=2048;
    ktsqr_multiplier=std::exp(std::log(maxktsqr/minktsqr)/N);
}

ChebyshevSolver::~ChebyshevSolver()
{
    delete[] wwork;
    delete[] ipwork;
    delete[] twork;

}
