/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
#include "chebyshev_amplitude.hpp"
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_math.h>

#include <cmath>
using std::abs;

const int MAXITER_VINT=1000;
const int MAXITER_UINT=MAXITER_VINT;
const REAL UVINTACCURACY=0.001;

enum BASIS_BOUNDARY_CONDITION
{
    CHEBYSHEV,
    CHEBYSHEV_ZERO1,        // Chebyshevs which are 0 at x=1
    CHEBYSHEV_ZERO2         // Chebyshevs which are 0 at x=\pm 1
};

const BASIS_BOUNDARY_CONDITION BC=CHEBYSHEV; //_ZERO2;



struct Integrand_helper
{
    ChebyshevAmplitudeSolver* N;
    REAL v;
    REAL y;
    unsigned int n;
    unsigned int m;
};

REAL Integrand_helperu(REAL v, void* p);
REAL Integrand_helperv(REAL v, void *p)
{
    Integrand_helper* h = (Integrand_helper*) p;
    h->v = v;
    
    gsl_function int_helper;
    int_helper.function=&Integrand_helperu;
    int_helper.params=h;
    REAL result, abserr;
    
    /*gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MAXITER_UINT);
    int status=gsl_integration_qag(&int_helper, -0.999, 0.999, 0, UVINTACCURACY, 
        MAXITER_UINT, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    */
    
    gsl_integration_cquad_workspace * workspace =gsl_integration_cquad_workspace_alloc(200);
    int status=gsl_integration_cquad(&int_helper,-0.999, 0.999, 0, UVINTACCURACY,
            workspace, &result, &abserr, NULL);
    gsl_integration_cquad_workspace_free(workspace);
    
    if (status) std::cerr << "Error " << status << " at " << LINEINFO
        << ": Result " << result << ", abserror: " << abserr 
        << " (v=" << v <<")" << endl;
    //cout << "u integral gave " << result << " at v=" << v << endl;
    return result;
    
    
}

REAL Integrand_helperu(REAL u, void* p)
{
    Integrand_helper* h = (Integrand_helper*) p;
    REAL result=0;
    unsigned int n=h->n;
    unsigned int m=h->m;
    REAL v=h->v;
    REAL m1 = h->N->M1(); REAL m2 = h->N->M2();
    if (std::abs(u-v)>0.00001)
    {
        result = h->N->Chebyshev(n, v) - std::exp(h->N->Delta(u,v))*h->N->Chebyshev(n, u);
        result /= std::abs(exp(h->N->Delta(u,v))-1);
    }
    else
        result = ( -(m1+m2)*(SQR(v)-1)*std::cos( n*std::acos(v) )
                + n*std::sqrt(1-SQR(v))*std::sin( n*std::acos(v) ) )
                / ( (m1+m2)*(SQR(v)-1) );
    
    result += h->N->Chebyshev(n, u)/std::sqrt(1.0+4.0*std::exp(-2.0*h->N->Delta(u,v)));
    result /= std::sqrt(1.0-SQR(u));
    result *= h->N->Chebyshev(m, u);

    return result;
    
}

void ChebyshevAmplitudeSolver::Solve(REAL maxy)
{
    Prepare();
    REAL alphabar = 0.2;

    /* Compute matrix F_{nm}
     * = \int_0^1 dv \int_-1^1 du (1-u^2)^(-1/2) T_m(u) {
     * (T_n(v) - exp(-delta)T_n(u) ) / Abs(exp(delta)-1)
     * + T_n(u)/Sqrt(1+4*exp(-2*delta)) }
     * here delta = (m1+m2)(u-v)
     */

    Integrand_helper h;
    h.N=this;
    h.y=0;

    gsl_function int_helper;
    int_helper.function=&Integrand_helperv;
    int_helper.params=&h;
    REAL result, abserr;

    std::vector< std::vector<REAL> > mat;

    for (unsigned int n=0; n<CHEBYSHEV_DEGREE; n++)
    {
        std::vector<REAL> tmpvec;
        #pragma omp parallel for
        for (unsigned int m=0; m<CHEBYSHEV_DEGREE; m++)
        {
            h.n=n;
            h.m=m;
    
            /*gsl_integration_workspace *workspace 
             = gsl_integration_workspace_alloc(MAXITER_UINT);
            int status=gsl_integration_qag(&int_helper, 0, 0.999, 0, UVINTACCURACY, 
                MAXITER_UINT, GSL_INTEG_GAUSS61, workspace, &result, &abserr);
            gsl_integration_workspace_free(workspace);
            */
            gsl_integration_cquad_workspace * workspace =gsl_integration_cquad_workspace_alloc(200);
            int status =gsl_integration_cquad(&int_helper, 0, 0.999, 0, UVINTACCURACY,
                workspace, &result, &abserr, NULL);
            gsl_integration_cquad_workspace_free(workspace);

            tmpvec.push_back(result);
            
            if (status) std::cerr << "Error " << status << " at " << LINEINFO
                << ": Result " << result << ", abserror: " << abserr << endl;

            cout << "#F_{" << n << "," << m << "}=" <<result << " relerr " << std::abs(abserr/result) << endl;
        }
        mat.push_back(tmpvec);
    }

    // Evolve up to maxy
    // Find maxyind corresponding to maxy
    unsigned int maxyind=YPoints();
    for (unsigned int i=1; i<=YPoints(); i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }


    // \pi_m \partial y a_m = (M_1+M_2) \alphabar \sum_n a_n f_{nm}
    for (unsigned int yind=1; yind<maxyind; yind++)
    {
        for (unsigned int aind=0; aind < mat.size(); aind++)
        {
            REAL newa=0;
            for (unsigned int tmpind=0; tmpind<mat.size(); tmpind++)
            {
                //cout << "coef[" << yind-1 << "][" << tmpind << "]*mat[tmpind][" << aind << endl;
                newa += coef[yind-1][tmpind]*mat[tmpind][aind];
            }
            if (aind==0) newa /= M_PI;
            else newa /= M_PI/2.0;
            newa *= (M1()+M2())*alphabar;
            newa *= yvals[yind]-yvals[yind-1];
            coef[yind][aind]=newa;

        }
    }

    ///DEBUG
    // Tulostetaan kun yind=5
    gsl_cheb_series *cs = gsl_cheb_alloc (CHEBYSHEV_DEGREE);
    for (unsigned int i=0; i<=CHEBYSHEV_DEGREE; i++)
    {        
        cs->c[i] = coef[5][i];
    }
    cs->a=0.0;
    cs->b=1.0;
    cs->order=CHEBYSHEV_DEGREE;

    for (unsigned int uind=0; uind<100; uind++)
    {
        REAL tmpu = uind/100.0;
        REAL ktsqr = Ktsqr(tmpu);
        cout << ktsqr << " " << gsl_cheb_eval(cs, tmpu) << endl;
    }

};

/*
 * Prepare everything
 * Allocate memory and compute the coefficients for the Chebyshev approx
 * of the initial condition
 */


REAL ICHelperf(REAL x, void* p)
{
    ICHelper* par = (ICHelper*) p;
    REAL ktsqr = par->N->Ktsqr(x);
    return par->N->InitialCondition(ktsqr);   
}

void ChebyshevAmplitudeSolver::Prepare()
{
    m1 = -std::log(MinKtsqr());
    m2 = std::log(MaxKtsqr());
    for (unsigned int yind=0; yind <= YPoints(); yind++)
    {
        REAL* tmpc = new REAL[CHEBYSHEV_DEGREE+1];
        coef.push_back(tmpc);
    }

    // Compute initial condition function expansion in the specified basis
    ICHelper help; help.N=this;
    for (unsigned int i=0; i<=CHEBYSHEV_DEGREE; i++)
    {
        // Evaluate coefficient a_i
        // = \int_-1^1 (1-x^2)^{-1/2} IC(x)*base[i](x)
        coef[0][i] = basis[i].DotProduct(ICHelperf, &help);
    }

    /*switch(BC)
    {
    case CHEBYSHEV:
    {

        gsl_cheb_series *cs = gsl_cheb_alloc (CHEBYSHEV_DEGREE);
        ICHelper help;
        help.N=this;
        gsl_function f;
        f.function = ICHelperf;
        f.params=&help;

        gsl_cheb_init(cs, &f, 0.0, 1.0);

        
        for (unsigned int i=0; i<=CHEBYSHEV_DEGREE; i++)
            coef[0][i]=cs->c[i];

        gsl_cheb_free(cs);
        break;
    }
    case CHEBYSHEV_ZERO1:
    case CHEBYSHEV_ZERO2:
        while(0){}
        
    }*/

    ComputeBasisVectors();

}

/*
 * Compute basis vectors
 * We set appropriate boundary conditions in order to make integrals easy
 */
void ChebyshevAmplitudeSolver::ComputeBasisVectors()
{
    // Make the basis orthonormal by Gram-Schmidt process

    switch(BC)
    {
    case CHEBYSHEV_ZERO2:
        {
        // Boundary condition: T_n(1)=0 and T_n(-1)=0
        ChebyshevVector v0(CHEBYSHEV_DEGREE);
        v0.SetComponent(2,1); v0.SetComponent(0,-1);    // T_2 -1
        v0.Normalize();
        basis.push_back(v0);
        cout << "Basis vector 0:" << v0 << endl;
        for (unsigned int j=3; j<=CHEBYSHEV_DEGREE; j++)
        {
            ChebyshevVector tmpvec(CHEBYSHEV_DEGREE);
            tmpvec.SetComponent(j, 1);
            if (GSL_IS_ODD(j))
                tmpvec.SetComponent(1,-1);
            else
                tmpvec.SetComponent(0,1);

            // Remove components in the directions of basis vectors
            for (unsigned int i=0; i<basis.size(); i++)
            {
                REAL dp = tmpvec.DotProduct(basis[i]);
                ChebyshevVector comp(CHEBYSHEV_DEGREE);
                comp=basis[i];
                comp = comp*dp;
                tmpvec = tmpvec - comp;
            }
            tmpvec.Normalize();
            cout << "Basis vector " << j-2 <<": " << tmpvec << endl;
            basis.push_back(tmpvec);
        }
        }
        break;
    case CHEBYSHEV:
        for (unsigned int i=0; i<=CHEBYSHEV_DEGREE; i++)
        {
            ChebyshevVector tmpvec(CHEBYSHEV_DEGREE);
            tmpvec.SetComponent(i, 1.0);
            tmpvec.Normalize();
            basis.push_back(tmpvec);
        }
        break;
    }
            


}


/*
 * Evaluate ith basis vector
 */
REAL ChebyshevAmplitudeSolver::Basis(unsigned int n, REAL x)
{
    if (n>=basis.size())
    {
        cerr<< "Asked basis vector "<< n << " but the largest index is " <<
        basis.size()-1 << " at " << LINEINFO << endl;
        return 0;
    }

    return basis[n].Evaluate(x);

}

ChebyshevAmplitudeSolver::~ChebyshevAmplitudeSolver()
{
    for (unsigned int yind=0; yind<YPoints(); yind++)
    {
        delete[] coef[yind];
    }
    gsl_cheb_free(cheb);

}

ChebyshevAmplitudeSolver::ChebyshevAmplitudeSolver()
{
    cheb = gsl_cheb_alloc (CHEBYSHEV_DEGREE);
    cheb->c = new REAL[CHEBYSHEV_DEGREE+1];
    cheb->a=-1.0; cheb->b=1.0;
    for (unsigned int i=1; i<=CHEBYSHEV_DEGREE; i++)
        cheb->c[i]=0;
    cheb->c[0]=2.0;
    oldn=0;
    
}

REAL ChebyshevAmplitudeSolver::M1()
{
    return m1;
}

REAL ChebyshevAmplitudeSolver::M2()
{
    return m2;
}

REAL ChebyshevAmplitudeSolver::Ktsqr(REAL u)
{
    // Comput ktsqr from u
    return std::exp( (M1()+M2())*u - M1() );
}

/*
 * Evaluate T_n(x)
 */
REAL ChebyshevAmplitudeSolver::Chebyshev(unsigned int n, REAL x)
{
    if (n > CHEBYSHEV_DEGREE) 
    {
        cerr << "Asked Chebyshev polynomial at degree n=" << n 
            << ", maximum  is " << CHEBYSHEV_DEGREE << endl;
        return -1;
    }
    if (oldn==n)
        return gsl_cheb_eval(cheb, x);
    cheb->c[oldn]=0.0;
    if (n==0)
    {
        cheb->c[0]=2.0;
    }
    else 
    {
        cheb->c[n]=1.0;
    }

    REAL result = gsl_cheb_eval_n(cheb, n, x);
    oldn=n;
    return result;
    
    
}

REAL ChebyshevAmplitudeSolver::Delta(REAL u, REAL v)
{
    return (M1()+M2())*(u-v);    
}
