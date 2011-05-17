/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "tools.hpp"
#include "chebyshev.hpp"
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_integration.h>
#include <cmath>

const unsigned int MAXITER_INNERPROD = 1000;
const REAL INNERPROD_ACCURACY = 0.0001;

/*
 * Compute the inner product between two Chebyshevs
 * E.g. integrate the product function over [-1,1] with weight
 * (1-x^2)^{-1/2}
 */
struct Inthelper_innerprod
{
    ChebyshevBasis* a;
    ChebyshevBasis* b;
    

};

REAL Inthelperf_innerprod(REAL x, void* p)
{
    Inthelper_innerprod *par = (Inthelper_innerprod*) p;
    return 1.0/std::sqrt(1.0-SQR(x))*par->a->Evaluate(x)*par->b->Evaluate(x);
}

REAL ChebyshevBasis::InnerProduct(ChebyshevBasis &vec)
{
    Inthelper_innerprod helper;
    helper.a=this;
    helper.b=&vec;

    gsl_function int_helper;
    int_helper.function=&Inthelperf_innerprod;
    int_helper.params=&helper;
    REAL result, abserr;
    /*
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MAXITER_INNERPROD);
    int status=gsl_integration_qag(&int_helper, -1, 1, 0,
        INNERPROD_ACCURACY, MAXITER_INNERPROD,
        GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    */
    gsl_integration_cquad_workspace * workspace =gsl_integration_cquad_workspace_alloc(200);
    int status =gsl_integration_cquad(&int_helper, -1, 1, 0, INNERPROD_ACCURACY,
         workspace, &result, &abserr, NULL);
    gsl_integration_cquad_workspace_free(workspace);

    return result;    


}

REAL ChebyshevBasis::Evaluate(REAL x)
{
    if (std::abs(x)>1)
    {
        cerr << "Asked the value of Chebyshev polynomial at x=" << x << " "
        << LINEINFO << endl;
        return 0;
    }

    return gsl_cheb_eval(cheb, x);

}

ChebyshevBasis::ChebyshevBasis(unsigned int d)
{
    degree = d;
    cheb = gsl_cheb_alloc (degree);
    cheb->a=-1.0;
    cheb->b=1.0;
    cheb->order=degree;
}

REAL ChebyshevBasis::Component(unsigned int n)
{
    if (n>degree)
    {
        cerr<< "Asked chebyshev-vector component " << n << ", but "
            << "dimension is " << degree << " at " << LINEINFO << endl;
        return 0.0;
    }

    return cheb->c[n];
}

void ChebyshevBasis::SetComponent(unsigned int c, REAL val)
{
    if (c>degree)
    {
        cerr<< "Tried to set component " << c << " to chebyshev vector but the "
            << "dimension is " << degree << " at " << LINEINFO << endl;
        return;
    }

    if (c==0) val*=2.0;
    cheb->c[c]=val;

}
