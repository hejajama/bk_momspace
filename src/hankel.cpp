/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

// Calculate discrete Hankel transformation using GSL
// So get N(r) from N(k) via
// N(r) = r^2 \int_0^\infty k*J_0(k*r)*N(k)

#include "hankel.hpp"
#include "config.hpp"
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <cmath>

struct Inthelper_hankel
{
    REAL r;
    Amplitude* N;
    REAL y;
};

/* Integrate \int d^2k/(2\pi) exp(-ik.r) N(k)
 * = \int dk k J_0(k*r) N(k)
 * = \int du exp(2*u) J_0[ exp(u)*r ]*N[ exp(u) ],
 * k=exp(u)
 */ 
REAL Inthelperf_hankel(REAL u, void* p)
{
    Inthelper_hankel* par = (Inthelper_hankel*)p;
    //REAL result = k*gsl_sf_bessel_J0(k* par->r) * par->N->N(SQR(k), par->y, false);
    REAL result = std::exp(2.0*u)*gsl_sf_bessel_J0(exp(u)*par->r)
        * par->N->N(std::exp(2.0*u), par->y, false);
    //cout << u << " " << std::exp(2.0*u) << " " << result << endl;
    return result;
}

void Hankel::PrintRAmplitude()
{
    cout << "# r  N" << endl;
/*
    for (int i=0; i<points; i++)
    {
        REAL tmpr = gsl_dht_k_sample(dht, i);
        cout << tmpr << " " << SQR(tmpr)*transformed[i] << endl;

    }
    return;
    */
    

    REAL minr = 5e-3; REAL maxr=1e3;
    REAL mult = std::pow(maxr/minr, 1.0/((REAL)points) );

    #pragma omp parallel for
    for (int i=0; i< points ; i++)
    {
        REAL r = minr*std::pow(mult, i);
        REAL result, abserr;
        
        Inthelper_hankel par;
        par.r=r; par.y=y, par.N=N;
        gsl_function int_helper;
        int_helper.function=&Inthelperf_hankel;
        int_helper.params=&par;

        const size_t maxiter = 300000;

        gsl_integration_workspace *workspace 
            = gsl_integration_workspace_alloc(maxiter);
        REAL minlnk = std::log(std::sqrt(N->Ktsqrval(0)));
        REAL maxlnk = std::log(std::sqrt(N->Ktsqrval(N->KtsqrPoints()-1)));
        int status = gsl_integration_qag(&int_helper, minlnk,
            maxlnk, 0, 0.02, maxiter, GSL_INTEG_GAUSS21, workspace,
            &result, &abserr);
        gsl_integration_workspace_free(workspace);

        if (status)
        {
            cerr << "Hankel transformation integral failed at " << LINEINFO <<
                ", r=" << r << ", result " << result << " relerr "
                << std::abs(abserr/result) << " y=" << y << endl;
        }
        
        cout << r << " " << SQR(r)*result << " " << result << endl;
    }
}


/*
 * Perform discrete Hankel transformation of the order 0
 */

Hankel::Hankel(Amplitude* amp, REAL y_, int npoints)
{
    y=y_;
    points=npoints;
    N=amp;
    
    sample = new REAL[points];
    transformed = new REAL[points];
/*
    REAL maxk = std::sqrt( N->KtsqrPoints()-2 );
    dht = gsl_dht_new(points, 0, maxk );
    
    for (int i=0; i<points; i++)
    {
        REAL tmpkt = gsl_dht_x_sample(dht, i);
        sample[i]=N->N(SQR(tmpkt), y);
    }

    gsl_dht_apply( dht, sample, transformed );
*/    
}

Hankel::~Hankel()
{
    delete[] sample;
    delete[] transformed;
    //gsl_dht_free(dht);
}
