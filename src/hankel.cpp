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
extern "C"
{
    #include "fourier/fourier.h"
}

struct Inthelper_hankel
{
    Amplitude* N;
    REAL y;
};

// Function to be transformed
REAL Helperf_hankel(REAL kt, void* p)
{
    Inthelper_hankel* par = (Inthelper_hankel*)p;
    //cout << "helper called at kt=" << kt << ", y=" << par->y << endl;
    return kt*par->N->N(SQR(kt), par->y);
}


/*
 * Perform Hankel transformation of the order 0
 */
REAL Hankel::Amplitude_r(REAL r, REAL y)
{
    Inthelper_hankel par;
    par.y=y; par.N=N;
    return SQR(r)*fourier_j0(r,Helperf_hankel,&par);
}

Hankel::Hankel(Amplitude* amp)
{
    N=amp;

    // Some initialisation stuff -- in principle nothing should be
    // changed here (1000 is the number of zeros of J0(x) that have been
    // encoded in the table in the file fourier.c)
    set_fpu_state();
    //gsl_set_error_handler_off();
    init_workspace_fourier(1000);
    set_fourier_precision(1.0e-12,1.0e-12);
       
}

Hankel::~Hankel()
{

}
