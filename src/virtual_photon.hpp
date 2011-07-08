#ifndef _VIRTUAL_PHOTON_H
#define _VIRTUAL_PHOTON_H

/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2011
 */

// Virtual photon wave function

#include "wave_function.hpp"
#include "config.hpp"

class VirtualPhoton : public WaveFunction
{
    public:
        VirtualPhoton();
        REAL PsiSqr_T_k(REAL Qsqr, REAL k, REAL z);
        REAL PsiSqr_L_k(REAL Qsqr, REAL k, REAL z);
        REAL Eps(REAL Qsqr, REAL z);

    private:
        // Quark masses and charges
        REAL mf[3];
        REAL ef[3];


};


#endif
