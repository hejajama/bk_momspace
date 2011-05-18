/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "tools.hpp"
#include "chebyshev.hpp"
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_integration.h>
#include <cmath>
#include <sstream>

const unsigned int MAXITER_INNERPROD = 1000;
const REAL INNERPROD_ACCURACY = 0.0001;

/*
 * Compute the inner product between two Chebyshevs/other functions
 * E.g. integrate the product function over [-1,1] with weight
 * (1-x^2)^{-1/2}
 */
 
REAL ChebyshevVector::DotProduct(ChebyshevVector &vec)
{
    if (Degree() != vec.Degree())
    {
        cerr << "Can't calclulate inner product between two vectors whose "
            << "dimensions are not the same! " << LINEINFO << endl;
        return 0;
    }

    REAL result=Component(0)*vec.Component(0)*M_PI;
    for (unsigned int n=1; n<=degree; n++)
        result += Component(n)*vec.Component(n)*M_PI/2.0;

    return result;

}


struct Inthelper_innerprod
{
    ChebyshevVector* a;
    REAL (*f)(REAL x, void* p) ;
    void* p;   

};

REAL Inthelperf_innerprod(REAL x, void* p)
{
    Inthelper_innerprod *par = (Inthelper_innerprod*) p;
    return 1.0/std::sqrt(1.0-SQR(x))*par->a->Evaluate(x)* (*(par->f))(x, par->p);
}

REAL ChebyshevVector::DotProduct(REAL(*f)(REAL x, void* p), void* p)
{

   Inthelper_innerprod helper;
    helper.a=this;
    helper.f=f;
    helper.p=p;

    gsl_function int_helper;
    int_helper.function=&Inthelperf_innerprod;
    int_helper.params=&helper;
    REAL result, abserr;

    gsl_integration_cquad_workspace * workspace =gsl_integration_cquad_workspace_alloc(200);
    int status =gsl_integration_cquad(&int_helper, -1, 1, 0, INNERPROD_ACCURACY,
         workspace, &result, &abserr, NULL);
    gsl_integration_cquad_workspace_free(workspace);

    return result;    


}

void ChebyshevVector::Normalize()
{
    REAL norm = std::sqrt(DotProduct((*this)));
    for (unsigned int i=0; i<=degree; i++)
    {
        SetComponent(i, Component(i)/norm );
    }
}

REAL ChebyshevVector::Evaluate(REAL x)
{
    if (std::abs(x)>1)
    {
        cerr << "Asked the value of Chebyshev polynomial at x=" << x << " "
        << LINEINFO << endl;
        return 0;
    }

    return gsl_cheb_eval(cheb, x);

}

ChebyshevVector::ChebyshevVector(unsigned int d)
{
    degree = d;
    cheb = gsl_cheb_alloc (degree);
    cheb->a=-1.0;
    cheb->b=1.0;
    cheb->order=degree;

    for (unsigned int i=0; i<=d; i++)
        SetComponent(i,0);
}

void ChebyshevVector::SetComponent(unsigned int c, REAL val)
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

REAL ChebyshevVector::Component(unsigned int c)
{
    if (c > degree)
    {
        cerr << "Asked component " << c << " which is larger than dimension "
            << degree << " " << LINEINFO << endl;
        return 0;
    }
    if (c==0) return cheb->c[c]/2.0;
    else return cheb->c[c];
}

unsigned int ChebyshevVector::Degree()
{
    return degree;
}

///// Operators
ChebyshevVector ChebyshevVector::operator+(ChebyshevVector &vec)
{
    if (Degree() != vec.Degree())
    {
        cerr << "Dimensio mismatch at " << LINEINFO << endl;
    }
    ChebyshevVector tmpvec(degree);
    for (unsigned int i=0; i<=degree; i++)
    {
        tmpvec.SetComponent(i, Component(i)+vec.Component(i));
    }
    return tmpvec;

}
ChebyshevVector ChebyshevVector::operator-(ChebyshevVector &vec)
{
    if (Degree() != vec.Degree())
    {
        cerr << "Dimensio mismatch at " << LINEINFO << endl;
    }
    ChebyshevVector tmpvec(degree);
    for (unsigned int i=0; i<=degree; i++)
    {
        tmpvec.SetComponent(i, Component(i)-vec.Component(i));
    }
    return tmpvec;

}
ChebyshevVector& ChebyshevVector::operator=(ChebyshevVector vec)
{
    if (Degree() != vec.Degree())
    {
        cerr << "Dimensio mismatch at " << LINEINFO << endl;
    }
    for (unsigned int i=0; i<=degree; i++)
    {
        SetComponent(i, vec.Component(i));
    }

    return *this;

}

ChebyshevVector& ChebyshevVector::operator*(REAL x)
{
    for (unsigned int i=0; i<=Degree(); i++)
    {
        cheb->c[i]*=x;
    }

    return *this;
}



std::ostream& operator<<(std::ostream& os, ChebyshevVector& v)
{
    std::stringstream ss;
    ss << "Chebyshev (";
    for (unsigned int i=0; i<v.Degree(); i++)
    {
        ss << v.Component(i) << ", ";
    }
    ss << v.Component(v.Degree()) << ")";
    return os << ss.str();

}
