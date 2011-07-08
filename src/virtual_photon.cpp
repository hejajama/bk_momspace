/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2011
 */

#include <cmath>
#include "virtual_photon.hpp"
#include "config.hpp"

REAL VirtualPhoton::PsiSqr_T_k(REAL Qsqr, REAL k, REAL z)
{
    return 0;
}


REAL VirtualPhoton::PsiSqr_L_k(REAL Qsqr, REAL k, REAL z)
{
    return 2.0*ALPHA_e/SQR(M_PI)*Qsqr*z*(1.0-z)*1.0;    ///TODO!
}

REAL VirtualPhoton::Eps(REAL Qsqr, REAL z)
{
    return 0;
}



VirtualPhoton::VirtualPhoton()
{
    ef[0]=2.0/3.0;      // u
    ef[1]=1.0/3.0;      // d
    ef[2]=1.0/3.0;      // s

}

