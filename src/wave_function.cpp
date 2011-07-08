/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "wave_function.hpp"
#include <iostream>

WaveFunction::WaveFunction()
{
    mode=VM_MODE_TOT;
}
/*
 * Integrate over z
 */

REAL WaveFunction::PsiSqr_T_intz_k(REAL Qsqr, REAL ksqr)
{
    return 0;    //TODO
}
REAL WaveFunction::PsiSqr_L_intz_k(REAL Qsqr, REAL ksqr)
{
    return 0;
}

/*
 * PsiSqr_intz
 * Returns |\Psi_L|^2, |\Psi_T|^2 or |\Psi_L|^2 + |\Psi_T|^2 depending
 * on the specified mode
 */
REAL WaveFunction::PsiSqr_intz_k(REAL Qsqr, REAL ksqr)
{
    switch (mode) 
    {
        case VM_MODE_TOT:
            std::cerr << "Can't calculate total wave function overlap in "
             << "WaveFunction::PsiSqr_intz!" << std::endl;
            break;
        case VM_MODE_L:
            return PsiSqr_L_intz_k(Qsqr, ksqr);
            break;
        case VM_MODE_T:
            return PsiSqr_T_intz_k(Qsqr, ksqr);
            break;
        default:
            std::cerr << "WaveFunctioN::PsiSqr_intz: Unknown mode " << mode 
                << std::endl;
     }
}  

REAL WaveFunction::PsiSqr_tot_k(REAL Qsqr, REAL ksqr, REAL z)
{
    return PsiSqr_T_k(Qsqr,ksqr,z)+PsiSqr_L_k(Qsqr,ksqr,z);
}

REAL WaveFunction::PsiSqr_tot_intz_k(REAL Qsqr, REAL ksqr)
{
    return PsiSqr_T_intz_k(Qsqr,ksqr)+PsiSqr_L_intz_k(Qsqr,ksqr);
}

void WaveFunction::SetMode(int m)
{
    mode=m;
}
