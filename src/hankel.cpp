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

REAL Inthelperf_hankel(REAL k, void* p)
{
    Inthelper_hankel* par = (Inthelper_hankel*)p;
    return k*gsl_sf_bessel_j0(k* par->r) * par->N->N(SQR(k), par->y, false);
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
    

    REAL minr = 1e-6; REAL maxr=1e1;
    REAL mult = std::pow(maxr/minr, 1.0/((REAL)points) );

    #pragma omp parallel for
    for (int i=0; i<points; i++)
    {
        REAL r = minr*std::pow(mult, i);
        REAL result, abserr;
        
        Inthelper_hankel par;
        par.r=r; par.y=y, par.N=N;
        gsl_function int_helper;
        int_helper.function=&Inthelperf_hankel;
        int_helper.params=&par;

        gsl_integration_cquad_workspace *workspace 
            = gsl_integration_cquad_workspace_alloc(1000);
        int status = gsl_integration_cquad(&int_helper, std::sqrt(N->Ktsqrval(0)),
            std::sqrt(N->Ktsqrval(N->KtsqrPoints()-2)), 0,
            0.02, workspace, &result, &abserr, NULL);
        gsl_integration_cquad_workspace_free(workspace);

        if (status)
        {
            cerr << "Hankel transformation integral failed at " << LINEINFO <<
                ", r=" << r << ", y=" << y << endl;
        }
        
        cout << r << " " << SQR(r)*result << endl;
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
