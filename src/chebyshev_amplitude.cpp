/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include "datafile.hpp"
#include "chebyshev_amplitude.hpp"
#include <vector>
#include <fstream>
#include <sstream>
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


const BASIS_BOUNDARY_CONDITION DEFAULT_BC=CHEBYSHEV;



struct Integrand_helper
{
    ChebyshevAmplitudeSolver* N;
    REAL u;
    REAL y;
    unsigned int n;
    unsigned int m;
};

REAL Integrand_helperv(REAL v, void* p);
REAL Integrand_helperu(REAL u, void *p)
{
    Integrand_helper* h = (Integrand_helper*) p;
    h->u = u;
    
    gsl_function int_helper;
    int_helper.function=&Integrand_helperv;
    int_helper.params=h;
    REAL result, abserr;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MAXITER_VINT);

    // Integrate over v with known "singular" point v=u
    // Actually it is not singularity there, only finite discontinuity
    // (at least in normal Chebyshev basis?)

    // Problems if u=\pm 1
    REAL maxy=1.0; REAL miny=-1.0;REAL pts[3];
    int status;
    if (u>1-1e-15 or u<-1+1e-15) {
        gsl_integration_cquad_workspace * workspace =gsl_integration_cquad_workspace_alloc(1000);
        status=gsl_integration_cquad(&int_helper,-1+1e-15, 1-1e-15, 0, UVINTACCURACY,
                workspace, &result, &abserr, NULL);
        gsl_integration_cquad_workspace_free(workspace);
    }
    else
    {
        pts[1]=u;
        pts[0]=miny; ; pts[2]=maxy;
        status = gsl_integration_qagp(&int_helper, pts, 3, 0, UVINTACCURACY,
            MAXITER_VINT, workspace, &result, &abserr);
    }
    gsl_integration_workspace_free(workspace); 
    /*int status=gsl_integration_qag(&int_helper, -0.999, 0.999, 0, UVINTACCURACY, 
        MAXITER_UINT, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    */
    
    /*gsl_integration_cquad_workspace * workspace =gsl_integration_cquad_workspace_alloc(200);
    int status=gsl_integration_cquad(&int_helper,-0.999, 0.999, 0, UVINTACCURACY,
            workspace, &result, &abserr, NULL);
    gsl_integration_cquad_workspace_free(workspace);
    */
    
    if (status) std::cerr << "Error " << status << " at " << LINEINFO
        << ": Result " << result << ", abserror: " << abserr 
        << " (u=" << u <<")" << endl;
    //cout << "u integral gave " << result << " at v=" << v << endl;
    return result;
    
    
}

REAL Integrand_helperv(REAL v, void* p)
{
    Integrand_helper* h = (Integrand_helper*) p;
    REAL result=0;
    unsigned int n=h->n;
    unsigned int m=h->m;
    REAL u=h->u;
     
    // We don't need to worry about the fact that the case u=v is not defined, as
    // the integration routine knows that something nasty is happening there
    if (std::abs(u-v)>1e-30)
    {
        result = h->N->Basis(n, v) - std::exp(h->N->Delta(u,v))*h->N->Basis(n, u);
        result /= std::abs(exp(h->N->Delta(u,v))-1);
    }
    else
    {
        cerr <<"u=" << u <<", v=" << v <<", first part of integral \approx zero! At "
         << LINEINFO << endl;
    }
    
    result += h->N->Basis(n, u)/std::sqrt(1.0+4.0*std::exp(-2.0*h->N->Delta(u,v)));
    //result /= std::sqrt(1.0-SQR(u));  <- This is inside qaws_table in u-integral
    result *= h->N->Basis(m, u);

    return result;
    
}


void ChebyshevAmplitudeSolver::SolveMatrix()
{
    /* Compute matrix F_{mn}
     * = \int_-1^1 dv \int_-1^1 du (1-u^2)^(-1/2) T_m(u) {
     * (T_n(v) - exp(-delta)T_n(u) ) / Abs(exp(delta)-1)
     * + T_n(u)/Sqrt(1+4*exp(-2*delta)) }
     * here delta = M(u-v)
     */



    

    mat.clear();
    for (unsigned int n=0; n<=chebyshev_degree; n++)
    {
        std::vector<REAL> tmpvec;
        for (unsigned int m=0; m<=chebyshev_degree; m++)
        {
            tmpvec.push_back(0.0);
        }
        mat.push_back(tmpvec);
    }
    

    int ready=0;
    #pragma omp parallel for
    for (unsigned int n=0; n<=chebyshev_degree; n++)
    {
        
        for (unsigned int m=0; m<=chebyshev_degree; m++)
        {
            REAL result, abserr;
            Integrand_helper h;
            h.N=this;
            h.y=0;
            gsl_function int_helper;
            int_helper.function=&Integrand_helperu;
            int_helper.params=&h;
            h.n=n;
            h.m=m;
    
            gsl_integration_workspace* workspace
            = gsl_integration_workspace_alloc(MAXITER_UINT);

            // Integrate with known singular ends \pm 1 caused by
            // weight (1-x^2)^(1/2) specified in qaws_table
            int status = gsl_integration_qaws( &int_helper, -1.0, 1.0, qaws_table,
                0, UVINTACCURACY, MAXITER_UINT, workspace, &result, &abserr);
            if (status )
            {
                cerr << "uint failed at " << LINEINFO << ", error " << status
                    << ", m=" << m <<", n=" << n <<", result: " << result << " relerr "
                    << std::abs(abserr/result) << endl;
            }

            mat[m][n]=result;
            ready++;

            if (ready%10==0)
                cout << "# ready " << ready << " / " << SQR(chebyshev_degree+1) << endl;
            
            if (status) std::cerr << "Error " << status << " at " << LINEINFO
                << ": Result " << result << ", relerror: "
                << std::abs(abserr/result) << endl;

            //cout << "#F_{m=" << m << ", n=" << n << "}=" <<result << " relerr " << std::abs(abserr/result) << endl;
        }


    }
    

}

/*
 * Save matrix into a file
 * Format is specified in README
 */
void ChebyshevAmplitudeSolver::SaveMatrix(std::string file)
{
    
    std::ofstream stream;
    stream.open(file.c_str());
    stream <<"# Solved with initial condition " << InitialConditionStr() << endl;
    stream <<"# Basis is " << boundary_condition << " (1=normal chebyshev, 2=boundary cond. E(1)=0"
        << endl;
    stream << "#M1=" << M1() << ", M2=" << M2() << endl;
    stream <<"###" << chebyshev_degree << endl;

    for (unsigned int m=0; m<=chebyshev_degree; m++)
    {
        stream << "#m=" << m << endl;
        for (unsigned int n=0; n<=chebyshev_degree; n++)
        {
            stream << mat[m][n] << endl;
        }
    }

    stream.close();
}

void ChebyshevAmplitudeSolver::LoadMatrix(std::string file)
{
    std::ifstream stream;
    stream.open(file.c_str());
    mat.clear();
    while (!stream.eof())
    {
        std::string line;
        std::getline(stream, line);

        // Dimension
        if (line.substr(0,3)=="###")
        {
            chebyshev_degree = StrToInt(line.substr(3,line.length()-3));
            break;
        }
    }

    for (unsigned int n=0; n<=chebyshev_degree; n++)
    {
        std::vector<REAL> tmpvec;
        for (unsigned int m=0; m<=chebyshev_degree; m++)
        {
            tmpvec.push_back(0.0);
        }
        mat.push_back(tmpvec);
    }

    unsigned int m=0;
    unsigned int n=0;

    // Coefficients are in order
    // F_{m=0,n=0} \n F_{m=0, n=1} \n ...

    while (!stream.eof())
    {
        std::string line;
        std::getline(stream, line);
        if (line.substr(0,1)=="#")
            continue;
        mat[m][n]=StrToReal(line);
        n++;
        if (n>chebyshev_degree)
        {
            m++; n=0;
        }
        if (m>chebyshev_degree)
            break;
    }
        
    if (m<chebyshev_degree)
        cerr << "We found only " << m << " different m values for the matrix" << endl;


}
 
REAL inthelperf_nonlin(REAL x, void* par);
struct inthelper_nonlin;

void ChebyshevAmplitudeSolver::Solve(REAL maxy)
{
    Prepare();
    REAL alphabar = 0.2;

    /* Compute matrix F_{nm}
     * = \int_-1^1 dv \int_-1^1 du (1-u^2)^(-1/2) T_m(u) {
     * (T_n(v) - exp(-delta)T_n(u) ) / Abs(exp(delta)-1)
     * + T_n(u)/Sqrt(1+4*exp(-2*delta)) }
     * here delta = M(u-v)
     */

    /*for ( int i=-100; i<100; i++)
    {
        REAL ktsqr = Ktsqr(i/100.0);
        cout << ktsqr << " "<< InitialCondition(ktsqr) << " " << N(ktsqr, 0) << endl;

    }
    return;
    */
    

    // Evolve up to maxy
    // Find maxyind corresponding to maxy
    unsigned int maxyind=YPoints();
    for (unsigned int i=1; i<=YPoints(); i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }


    // \partial y a_m = (M_1+M_2) \alphabar \sum_n a_n f_{nm}
    for (unsigned int yind=1; yind<maxyind; yind++)
    {
        for (unsigned int aind=0; aind < mat.size(); aind++)
        {
            REAL dera=0;
            for (unsigned int tmpind=0; tmpind<mat.size(); tmpind++)
            {
                //cout << "coef[" << yind-1 << "][" << tmpind << "]*mat[tmpind][" << aind << endl;
                dera += coef[yind-1][tmpind]*mat[aind][tmpind];
            }
            dera *= (M1()+M2())*alphabar;
            dera -= alphabar*NonLinear(aind, yind-1);
            coef[yind][aind]=coef[yind-1][aind] + (yvals[yind]-yvals[yind-1])*dera;

        }
    }

    ///DEBUG
    // Tulostetaan kun yind=5
    gsl_cheb_series *cs = gsl_cheb_alloc (chebyshev_degree);
    for (unsigned int i=0; i<=chebyshev_degree; i++)
    {        
        cs->c[i] = coef[0][i];
    }
    cs->a=0.0;
    cs->b=1.0;
    cs->order=chebyshev_degree;

    for (int uind=-100; uind<100; uind++)
    {
        REAL tmpu = uind/100.0;
        REAL ktsqr = Ktsqr(tmpu);
        REAL amp=0; REAL amp2=0; REAL amp3=0;
        for (unsigned int d=0; d<=chebyshev_degree; d++)
        {
            amp += coef[0][d]*basis[d].Evaluate(tmpu);
            amp2 += coef[1][d]*basis[d].Evaluate(tmpu);
            amp3 += coef[5][d]*basis[d].Evaluate(tmpu);
        }
        
        cout << ktsqr << " " << amp << " " << amp2 << " " << amp3 << endl;
    }

};

/*
 * Compute \int du/Sqrt[1-u^2] E_m(x) N(x)^2
 */

/*
 * Contribution from the non-linear term -N(x)^2
 * Can be computed using the properties of the Chebyshev polynomials
 * m is the index of the coefficient which evolution we are computing,
 * meaning that integratin \int dx (1-x^2)^{-1/2} E_m(x) (\sum_i a_i E_i(x) )^2
 */
REAL Pi(unsigned int i)
{
    if (i==0) return M_PI;
    else return M_PI/2.0;
}

struct inthelper_nonlin
{
    ChebyshevAmplitudeSolver* N;
    unsigned int m;
    REAL y;
};

REAL inthelperf_nonlin(REAL x, void* par)
{
    inthelper_nonlin* p = (inthelper_nonlin*) par;
    return p->N->Basis(p->m, x)*SQR(p->N->N(p->N->Ktsqr(x), p->y) );
    
    
}

REAL ChebyshevAmplitudeSolver::NonLinear(unsigned int m, unsigned int yind)
{
    REAL result=0;

    // Weight function (1-x^2)^(-1/2) = [ (1-x)*(1+x) ]^(-1/2)
    // -> (x-a)^\alpha (b-x)^\beta \log^\mu (x-a) \log^\nu (b-x) weight function,
    // a=-1, b=1 => \alpha=-1/2, \beta=-1/2, \mu=0, \nu=0
    gsl_integration_qaws_table* table; table = gsl_integration_qaws_table_alloc( -0.5, -0.5, 0, 0);
    const int MAXITER_NONLINEAR = 1000;
    const REAL NONLINEAR_ACCURACY = 0.001;
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(MAXITER_NONLINEAR);
    inthelper_nonlin helper;
    helper.N=this;
    helper.m=m; helper.y=yvals[yind];
    gsl_function int_helper;
    int_helper.function=inthelperf_nonlin;
    int_helper.params=&helper;
    REAL abserr;
    int status = gsl_integration_qaws( &int_helper, -1.0, 1.0, table,
        0, NONLINEAR_ACCURACY, MAXITER_NONLINEAR, workspace, &result, &abserr);
    if (status)
    {
        cerr << "Error " << status << " at " << LINEINFO << ": result " << result
            << ", relerr " << std::abs(abserr/result) << ", m=" << m << ", yind=" << yind
            << endl;
    }
    gsl_integration_workspace_free(workspace);
    gsl_integration_qaws_table_free(table);
    return result;
    
    switch(boundary_condition)
    {
        case CHEBYSHEV:


            break;
        case CHEBYSHEV_ZERO:
            
            /*
            for (unsigned int n=0; n<=chebyshev_degree; n++)
            {
                for (unsigned int nn=0; nn<=chebyshev_degree; nn++)
                {
                    REAL tmpres=0;
                    for (unsigned int i2=0; i2<=chebyshev_degree; i2++)
                    {
                        for (unsigned int i3=0; i3<=chebyshev_degree; i3++)
                        {
                            if (i2+i3 <=chebyshev_degree)
                                tmpres += BasisVector(m).Component(i2+i3)
                                    * Pi(i2+i3);
                            
                            int abs = std::abs(static_cast<int>(i2)-static_cast<int>(i3));
                            tmpres += BasisVector(m).Component(abs)
                                * Pi(abs);
                            tmpres *= BasisVector(n).Component(i2)
                                    * BasisVector(nn).Component(i3);
                        }
                    }
                    tmpres *= coef[yind][n]*coef[yind][nn];
                    tmpres *= 0.5;
                    result +=tmpres;
                }
            }
            */
                            

            break;
    }
    return result;

}

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
extern "C"{
void dfct(int, double *, double *, int *, double *);
};
void ChebyshevAmplitudeSolver::Prepare()
{
    m1 = -std::log(MinKtsqr());
    m2 = std::log(MaxKtsqr());
    m1=0;   // => v \in [-1,1]

    
    for (unsigned int yind=0; yind <= YPoints(); yind++)
    {
        std::vector<REAL> tmpvec;
        for (unsigned int i=0; i<=chebyshev_degree; i++)
        {
            tmpvec.push_back(0.0);
        }
        coef.push_back(tmpvec);
    }
    ComputeBasisVectors();

    // Compute initial condition function expansion in the specified basis
    ICHelper help; help.N=this;
    coef[0][0]=0;
    for (unsigned int i=0; i<=chebyshev_degree; i++)
    {
        // Evaluate coefficient a_i
        // = \int_-1^1 (1-x^2)^{-1/2} IC(x)*base[i](x)
        coef[0][i] = basis[i].DotProduct(ICHelperf, &help);
    }

}

/*
 * Compute basis vectors
 * We set appropriate boundary conditions in order to make integrals easy
 */
void ChebyshevAmplitudeSolver::ComputeBasisVectors()
{
    // Make the basis orthonormal by Gram-Schmidt process
    switch(boundary_condition)
    {        
    case CHEBYSHEV_ZERO:
    {
        // Boundary condition: T_n(1)=0
        basis.clear();
        ChebyshevVector v1(chebyshev_degree);
        v1.SetComponent(1,1.0); v1.SetComponent(0,-1.0);    // T_1 -1
        v1.Normalize();
        basis.push_back(v1);
        for (unsigned int j=2; j<=chebyshev_degree; j++)
        {
            ChebyshevVector tmpvec(chebyshev_degree);
            tmpvec.SetComponent(j, 1.0);
            tmpvec.SetComponent(0, -1.0);
            // Remove components in the directions of basis vectors
            ChebyshevVector comp(chebyshev_degree);
            for (unsigned int i=0; i<basis.size(); i++)
            {
                REAL dp = tmpvec.DotProduct( BasisVector(i) );
                comp=basis[i];
                comp = comp*dp;
                tmpvec = tmpvec - comp;
            }
            tmpvec.Normalize();
            //cout << "Basis vector " << j-2 <<": " << tmpvec << endl;
            basis.push_back(tmpvec);
        }
        // Let's allways have chebyshev_degree+1 vectors in the basis vector
        ChebyshevVector v0(chebyshev_degree);
        basis.push_back(v0);
        break;
    }
    case CHEBYSHEV:
        for (unsigned int i=0; i<=chebyshev_degree; i++)
        {
            ChebyshevVector tmpvec(chebyshev_degree);
            tmpvec.SetComponent(i, 1.0);
            tmpvec.Normalize();
            basis.push_back(tmpvec);
        }
        break;
    }
            
}

/*
 * N
 * Evaluate the scattering amplitude at given ktsqr and y
 * So evaluates the amplitude in the basis used
 * Interpolates between rapiditeis somehow (TODO)
 */
REAL ChebyshevAmplitudeSolver::N(REAL ktsqr, REAL y)
{
    REAL yind=0;    ///TODO
    REAL u = U(ktsqr);

    // Efective coefficients (interpolate with spline?)
    REAL cf[chebyshev_degree+1];
    for (unsigned int i=0; i<=chebyshev_degree; i++)
    {
        cf[i] = coef[0][i];
    }

    REAL result=0;
    for (unsigned int i=0; i<=chebyshev_degree; i++)
    {
        result += cf[i]*Basis(i, u);
    }
    return result;
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

    gsl_cheb_free(cheb);

}

ChebyshevAmplitudeSolver::ChebyshevAmplitudeSolver()
{
    chebyshev_degree=DEFAULT_CHEBYSHEV_DEGREE;
    boundary_condition=DEFAULT_BC;
    cheb = gsl_cheb_alloc (chebyshev_degree);
    cheb->c = new REAL[chebyshev_degree+1];
    cheb->a=-1.0; cheb->b=1.0;
    for (unsigned int i=1; i<=chebyshev_degree; i++)
        cheb->c[i]=0;
    cheb->c[0]=2.0;
    oldn=0;

    // Weight function for dot products (1-x^2)^(-1/2) = [ (1-x)*(1+x) ]^(-1/2)
    // -> (x-a)^\alpha (b-x)^\beta \log^\mu (x-a) \log^\nu (b-x) weight function,
    // a=-1, b=1 => \alpha=-1/2, \beta=-1/2, \mu=0, \nu=0
    qaws_table = gsl_integration_qaws_table_alloc( -0.5, -0.5, 0, 0);
    
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

REAL ChebyshevAmplitudeSolver::U(REAL ktsqr)
{
    // u from ktsqr
    return M1() + std::log(ktsqr)/( M1()+M2() );
}
/*
 * Evaluate T_n(x)
 */
REAL ChebyshevAmplitudeSolver::Chebyshev(unsigned int n, REAL x)
{
    cerr << "ChebyshevAmplitudeSolver::Chebyshev called!" << endl;
    if (n > chebyshev_degree) 
    {
        cerr << "Asked Chebyshev polynomial at degree n=" << n 
            << ", maximum  is " << chebyshev_degree << endl;
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

void ChebyshevAmplitudeSolver::SetChebyshevDegree(unsigned int d)
{
    chebyshev_degree=d;
}

ChebyshevVector& ChebyshevAmplitudeSolver::BasisVector(unsigned int n)
{
    if (n>=basis.size())
    {
        cerr << "Asked basis vector " << n << " but maximum index is "
        << basis.size()-1 << endl;
        ChebyshevVector *tmp = new ChebyshevVector(chebyshev_degree);
        return (*tmp);
    }
    return basis[n];
}

/*
 * Sets boundary condition
 * Notice: Doesn't calculate new basis vectors!
 * So method ComputeBasisVectors() should be called separately
 */
void ChebyshevAmplitudeSolver::SetBoundaryCondition(BASIS_BOUNDARY_CONDITION bc)
{
    boundary_condition=bc;
}

REAL ChebyshevAmplitudeSolver::MinKtsqr()
{
    return std::exp(-std::log(maxktsqr));
}


REAL ChebyshevAmplitudeSolver::MaxKtsqr()
{
    return maxktsqr;
}
