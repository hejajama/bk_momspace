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
#include <iomanip>
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

const int MAXITER_KIN=1000;
const REAL KINACCURACY=0.01;


const BASIS_BOUNDARY_CONDITION DEFAULT_BC=CHEBYSHEV_ZERO;



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

    REAL maxy=1.0; REAL miny=-1.0;REAL pts[3];
    int status;

    if (h->N->KinematicalConstraint()==true)
    {
        // Integration limits are [-1,u]

        if (u < -1.0+1e-15)
            return 0;
        
        gsl_integration_cquad_workspace *workspace 
        = gsl_integration_cquad_workspace_alloc(MAXITER_VINT);
        status = gsl_integration_cquad(&int_helper, -1.0, u-1e-15, 0,
            UVINTACCURACY, workspace, &result, &abserr, NULL);
        gsl_integration_cquad_workspace_free(workspace);
    }

    else    // No kinematical constraint, integrate over v \in [-1,1]
    {
        // Integrate over v with known "singular" point v=u
        // Actually it is not singularity there, only finite discontinuity
        // (at least in normal Chebyshev basis?)

        // Problems if u=\pm 1
        if (u>1-1e-15 or u<-1+1e-15 ) {
            gsl_integration_cquad_workspace * workspace =gsl_integration_cquad_workspace_alloc(1000);
            status=gsl_integration_cquad(&int_helper,-1+1e-15, 1-1e-15, 0, UVINTACCURACY,
                    workspace, &result, &abserr, NULL);
            gsl_integration_cquad_workspace_free(workspace);
        }
        else
        {
            gsl_integration_workspace *workspace 
            = gsl_integration_workspace_alloc(MAXITER_VINT);
            pts[1]=u;
            pts[0]=miny; ; pts[2]=maxy;
            status = gsl_integration_qagp(&int_helper, pts, 3, 0, UVINTACCURACY,
                MAXITER_VINT, workspace, &result, &abserr);
            gsl_integration_workspace_free(workspace);
        }
     }
    
    if (status) std::cerr << "Error " << status << " at " << LINEINFO
        << ": Result " << result << ", abserror: " << abserr 
        << " (u=" << u <<")" << endl;
    //cout << "u integral gave " << result << " at v=" << v << endl;
    if (h->N->RunningCoupling()==PARENT_DIPOLE)
        result *= Alphabar_s(h->N->Ktsqr(u));     ///TODO: ok?
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
     *
     * If kinematical constraint is applied, limits for vint are
     * [-1,u], case v>0 is coevered by Kinematic() function
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
    for (int m=0; m<=chebyshev_degree; m++)
    {
        
        for (unsigned int n=0; n<=chebyshev_degree; n++)
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
    stream << std::setprecision(9) << "#M1=" << M1() << ", M2=" <<
        std::setprecision(9) << M2() << ", minktsqr=" << MinKtsqr()
            << ", maxktsqr=" << MaxKtsqr() << endl;
    stream <<"###" << chebyshev_degree << endl;
    stream <<"###" << std::setprecision(9) << M2() << endl;

    for (unsigned int m=0; m<=chebyshev_degree; m++)
    {
        stream << "#m="  << m << endl;
        for (unsigned int n=0; n<=chebyshev_degree; n++)
        {
            stream <<  std::scientific << std::setprecision(15) << mat[m][n] << endl;
        }
    }

    stream.close();
}

void ChebyshevAmplitudeSolver::LoadMatrix(std::string file)
{
    std::ifstream stream;
    stream.open(file.c_str());
    mat.clear();
    int setting=0;
    while (!stream.eof())
    {
        std::string line;
        std::getline(stream, line);

        // Dimension
        if (line.substr(0,3)=="###")
        {
            if (setting==0)
                chebyshev_degree = StrToInt(line.substr(3,line.length()-3));
            if (setting==1)
            {
                m2 = StrToReal(line.substr(3, line.length()-3));
                maxktsqr = std::exp(m2);
                minktsqr = std::exp(-m2);
            }
            setting++;
            if (setting>1)
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
    REAL alphabar = 0.2;

    /* Evolution is derived by matrix F_{nm}
     * = \int_-1^1 dv \int_-1^1 du (1-u^2)^(-1/2) T_m(u) {
     * (T_n(v) - exp(-delta)T_n(u) ) / Abs(exp(delta)-1)
     * + T_n(u)/Sqrt(1+4*exp(-2*delta)) }
     * here delta = M(u-v)
     * which is computed beforehand. If kinematical constraint is applied,
     * this integral contains extra constraint \theta(u-v)
     */

    if (YPoints()==0)
    {
        int ypoints = static_cast<int>(maxy/delta_y) + 1;
        yvals.clear();
        for (int i=0; i<=ypoints; i++)
        {
            yvals.push_back(i*delta_y);
        }
        Prepare();
    }

    // Evolve up to maxy
    // Find maxyind corresponding to maxy
    unsigned int maxyind=YPoints();
    for (unsigned int i=1; i<=YPoints(); i++)
    {
        if (yvals[i]>maxy)
            { maxyind=i; break; }
    }
    cout << "# Evolving up to y=" << yvals[maxyind] << ", M=" << M() << endl;

    // \partial y a_m = (M_1+M_2) \alphabar \sum_n a_n f_{nm}
    for (unsigned int yind=1; yind<=maxyind; yind++)
    {
        cerr << "# y=" << yvals[yind] << endl;
        #pragma omp parallel for
        for (int aind=0; aind <= chebyshev_degree; aind++)
        {
            REAL dera=0;
            for (unsigned int tmpind=0; tmpind<= chebyshev_degree; tmpind++)
            {
                //cout << "coef[" << yind-1 << "][" << tmpind << "]*mat[tmpind][" << aind << endl;
                dera += coef[yind-1][tmpind]*mat[aind][tmpind];
            }
            if (RunningCoupling()==CONSTANT)
            {
                dera *= (M1()+M2())*alphabar;
                dera -= alphabar*NonLinear(aind, yind-1);
            }
            else
            {
                dera *= (M1()+M2());
                dera -= NonLinear(aind, yind-1);
            }

            if (kinematic_constraint)
            {
                REAL kc = Kinematic(aind, yind-1);
                dera += alphabar*kc;    
                /*** LO
                for (unsigned int tmpind=0; tmpind <= chebyshev_degree; tmpind++)
                {
                    kc += Kinematic(aind, tmpind, yind-1);
                }
                kc *= alphabar;
                dera -= kc;
                */
            }

            REAL newcoef = coef[yind-1][aind] + (yvals[yind]-yvals[yind-1])*dera;
            if (std::abs( (coef[yind-1][aind] - newcoef) / coef[yind-1][aind]) > 0.1)
            {
                cerr << "Large change from " << coef[yind-1][aind] << " to " << newcoef
                    << " at aind=" << aind <<", yind=" << yind << endl;
            }
            
            coef[yind][aind]=coef[yind-1][aind] + (yvals[yind]-yvals[yind-1])*dera;

        }
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
    REAL result = p->N->Basis(p->m, x)*SQR(p->N->N(p->N->Ktsqr(x), p->y) );
    if (p->N->RunningCoupling()==PARENT_DIPOLE)
        result *= Alphabar_s(p->N->Ktsqr(x));
    return result; 
    
    
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
 * Contribution from the kinematical constraint in lowest order
 * Assume that a(Y) is a smooth, "slowly" varying function
 * => a(Y - M(v-u)) = a(Y) - M(v-u) a'(Y)
 * a'(Y) is calculated using linear approximation
 * we assume that a(x)=0 for all x<0
 *
 * Thus we integrate
 * \int_{-1}^1 du (1-u^2)^{-1/2} T_m(u) \int_u^1 dv M(v-u)*T_n(v) / Abs[1-exp[M(u-v)]]
 */

struct Inthelper_kinematic
{
    unsigned int n,m;
    REAL u;
    ChebyshevAmplitudeSolver* N;
    REAL y;
    uint yind;
    REAL sum_anTn;  // \sum_n a_n(Y) T_n(u)
};

REAL Inthelperf_kinematic_v_lo(REAL v, void* p)
{
    Inthelper_kinematic* par = (Inthelper_kinematic*)p;
    REAL result=0;
    
    if (std::abs(par->u - v) < 1e-10)
        result = 1.0;
    else
        result = (v - par->u)/std::abs( 1.0-std::exp( par->N->M2()*(par->u-v) ) );
    result *= par->N->Basis(par->n, v);
    return result;

}


/*
 * Contribution from the kinematical constraint exactly
 * In this case the rapiditiy-evolution matrix is computed with
 * additional constraint \theta(u-v)
 *
 * This function computes the contribution which can't be calculated
 * only once and saved in one matrix
 * Namely integrates
 * \int_{-1}^1 du/Sqrt[1-u^2] T_m(u) \int_u^1 dv \sum_n [
 *  (a_n(Y-M(u-v))T_n(v) - exp(M(u-v))*T_n(u)*a_n(Y) ) /
 *  |1-exp(M(u-v))| + T_n(u)/Sqrt[1+4*exp(-2*M*(u-v)) ] ]
 *
 * 1/Sqrt[1-u^2] has been taken into account in integral over u
 */



REAL Inthelperf_kinematic_v(REAL v, void* p)
{
    Inthelper_kinematic* par = (Inthelper_kinematic*)p;
    REAL result=0;

    if (std::abs(par->u - v) < 1e-15)
        return 0.0;

    // \sum_n a_n(Y-M(v-u)) T_n(v)
    // Note: Coefficient(y<0) = Coefficient(y=0)
    REAL nonlocalsum=0;
    for (uint i=0; i<=par->N->ChebyshevDegree(); i++)
    {
        nonlocalsum += par->N->CoefficientY(par->y - par->N->M()*(v - par->u), i )
                    * par->N->Basis(i, v);
    }

    /*if (std::abs(par->u-v)<1e-6)
        cerr << "u=" << par->u << ", v=" << v <<", sum=" << par->sum_anTn
            <<", nonlocsum=" << nonlocalsum << ", y=" << par->y << endl;
*/
    result = (nonlocalsum - std::exp( par->N->M()*(par->u - v) )*par->sum_anTn)
        /std::abs( 1.0-std::exp(par->N->M()*(par->u - v)) )
        + par->sum_anTn/std::sqrt( 1.0 + 4.0*std::exp(-2.0*par->N->M()*(par->u - v)) );
    return result;


}

REAL Inthelperf_kinematic_u(REAL u, void * p)
{
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(MAXITER_KIN);
    Inthelper_kinematic* par = (Inthelper_kinematic*)p;
    par->u=u;
    gsl_function int_helper;
    int_helper.function=&Inthelperf_kinematic_v;
    int_helper.params=p;

    //// Not needed in lo
    if (u> 1.0 - 1e-15) return 0;
    REAL sum = 0;
    for (uint i=0; i<=par->N->ChebyshevDegree(); i++)
    {
        sum += par->N->Coefficient(par->yind, i) * par->N->Basis(i, u);
    }
    par->sum_anTn=sum;
    ////

    REAL result, abserr;

    int status = gsl_integration_qag(&int_helper, u+1e-15, 1.0, 0, KINACCURACY,
        MAXITER_KIN, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

    if (status)
    {
        cerr << "v-integration failed at " << LINEINFO << ", code " << status
            << ", u=" << u /* << ", n=" << ((Inthelper_kinematic*)p)->n */
            << ", m=" << ((Inthelper_kinematic*)p)->m << ", result: " << result
            << ", relerr " << std::abs(abserr/result) << endl;
            /*
            for (REAL x=u; x<=1.0; x+=0.001)
                cout << x << " " << Inthelperf_kinematic_v(x, par) << endl;
            exit(1);
            */
    }
    
    gsl_integration_workspace_free(workspace);

    result *= ((Inthelper_kinematic*)p)->N->Basis(((Inthelper_kinematic*)p)->m, u);
    
    return result;
}
 
REAL ChebyshevAmplitudeSolver::Kinematic(unsigned int m,
    unsigned int yind)
{
    
    REAL result, abserr;
    // Weight function (1-x^2)^(-1/2) = [ (1-x)*(1+x) ]^(-1/2)
    // -> (x-a)^\alpha (b-x)^\beta \log^\mu (x-a) \log^\nu (b-x) weight function,
    // a=-1, b=1 => \alpha=-1/2, \beta=-1/2, \mu=0, \nu=0
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(MAXITER_KIN);
    gsl_integration_qaws_table *table = gsl_integration_qaws_table_alloc( -0.5, -0.5, 0, 0);
    
    Inthelper_kinematic helper;
    helper.N=this;
    helper.m=m;
    helper.y=yvals[yind];
    helper.yind=yind;

    gsl_function int_helper;
    int_helper.function=&Inthelperf_kinematic_u;
    int_helper.params=&helper;

    int status = gsl_integration_qaws( &int_helper, -1.0, 1.0, table,
        0, KINACCURACY, MAXITER_KIN, workspace, &result, &abserr);

    if (status and result>1e-3 )
    {
        cerr << "u-integration failed at " << LINEINFO << ", code " << status
            /* << ", n=" << n */
            << ", m=" << m << ", result: " << result
            << ", relerr " << std::abs(abserr/result) << endl;
    }

    gsl_integration_workspace_free(workspace);
    gsl_integration_qaws_table_free(table);

    result *= M();

    //// LO ONLY
    // Then multiply by a_n'(Y)
    // REAL dera = (coef[yind][n] - coef[yind-1][n])/(yvals[yind]-yvals[yind-1]);
    // result *= dera;
    
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
    cout << "# IC coef. (";
    for (unsigned int i=0; i<=chebyshev_degree; i++)
    {
        // Evaluate coefficient a_i
        // = \int_-1^1 (1-x^2)^{-1/2} IC(x)*base[i](x)
        coef[0][i] = basis[i].DotProduct(ICHelperf, &help);
        cout << coef[0][i] << ", " ;
    }
    cout << ")" << endl;


    ///DEBUG
   /*for (REAL u=-1; u<=1; u+=0.01)
    {
        REAL tmpres=0;
        for (int i=0; i<=chebyshev_degree; i++)
            tmpres += coef[0][i]*Basis(i, u);
        cout << Ktsqr(u) << " " << tmpres << " " << InitialCondition(Ktsqr(u)) << endl; 
    }
    exit(1); */
    

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
 * Interpolates between rapidites linearly (TODO: spline?)
 * Or interpolate coefficients?
 */
REAL ChebyshevAmplitudeSolver::N(REAL ktsqr, REAL y, bool bspline, bool derivative)
{
    // Find yind
    // val[index]  is smaller than (or equal) y
    int yind=-1;
    for (unsigned int i=0; i<yvals.size()-1; i++)
    {
        if (yvals[i]<=y and yvals[i+1]>y)
        {
            yind=i;
            break;
        }
    }
    if (yind<0)
        yind=yvals.size()-1;
        
    REAL u = U(ktsqr);

    // Interpolate linearly

    REAL cf0[chebyshev_degree+1];   // at lower rapidity
    REAL cf1[chebyshev_degree+1];   // at higher rapidity
    for (unsigned int i=0; i<=chebyshev_degree; i++)
    {        
        // interpolate coefficients:
        cf0[i] = coef[yind][i];
        if (yind < yvals.size()-1)
            cf1[i] = coef[yind+1][i];   
    }

    REAL n0=0; REAL n1=0;
    for (unsigned int i=0; i<=chebyshev_degree; i++)
    {
        n0 += cf0[i]*Basis(i, u);
        n1 += cf1[i]*Basis(i,u);
    }

    // Interpolate
    if (y - yvals[yind] < 0.00001)  // Don't interpolate
        return n0;
    if (yind >= yvals.size()-1) // Extrapolate
    {
        cerr << "Should extrapolate, as y=" << y<<" , falling back to y=" <<
        yvals[yvals.size()-1] << endl;
        return n1;
    }
    else
    {
        return n0 + (y-yvals[yind])*(n1-n0)/(yvals[yind+1]-yvals[yind]);
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

REAL ChebyshevAmplitudeSolver::M()
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

unsigned int ChebyshevAmplitudeSolver::ChebyshevDegree()
{
    return chebyshev_degree;
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

REAL ChebyshevAmplitudeSolver::Coefficient( unsigned int yind, unsigned int degree)
{
    return coef[yind][degree];
}


/*
 * Coefficient a(Y) as a function of Y
 * Interpolates linearly in Y
 * if Y<0 returns a(0)
 */

REAL ChebyshevAmplitudeSolver::CoefficientY(REAL y, uint degree)
{
    if (y<0.00001)
        return Coefficient(0, degree);
    
    // Find yind so that yvals[yind]<y and yvals[yind+1]>y
    int yind=-1;
    for (uint i=0; i<yvals.size()-1; i++)
    {
        if (yvals[i]<=y and yvals[i+1]>y)
        {
            yind=i;
            break;
        }
    }
    if (yind==-1)   // Didn't find, fall back to y=yind[yvals-1]
    {
        if (std::abs(y - yvals[yvals.size()-1]) > 0.05)
        {
            cerr << LINEINFO  << ": Asked Coefficient of degree " << degree << " at too large"
                << " value of y=" << y <<", falling back to y="
                << yvals[yvals.size()-1] << endl;
        }
        return Coefficient(yvals.size()-1, degree);

    }

    // Interpolate linearly
    return Coefficient(yind, degree) +
        (y-yvals[yind])* ( Coefficient(yind+1, degree) - Coefficient(yind, degree) )
        / (yvals[yind+1] - yvals[yind] ) ;
    

}
