#ifndef _WAVE_FUNCTION_H
#define _WAVE_FUNCTION_H

/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2011
 */

#include <string>
#include "config.hpp"

// General abstraction class for wave functions

const int VM_MODE_TOT=1; const int VM_MODE_L=2; const int VM_MODE_T=3;

// _k postfix <=> k space wavefunction, similarly _r for coord. space

class WaveFunction{
    public:
        WaveFunction();
        virtual REAL PsiSqr_T_k(REAL Qsqr, REAL ksqr, REAL z) = 0;
        virtual REAL PsiSqr_L_k(REAL Qsqr, REAL ksqr, REAL z) = 0;
        REAL PsiSqr_T_intz_k(REAL Qsqr, REAL ksqr);
        REAL PsiSqr_L_intz_k(REAL Qsqr, REAL ksqr);
        //virtual std::string GetParamString();
        REAL PsiSqr_tot_k(REAL Qsqr, REAL ksqr, REAL z);
        REAL PsiSqr_tot_intz_k(REAL Qsqr, REAL ksqr);
        REAL PsiSqr_intz_k(REAL Qsqr, REAL ksqr);
        void SetMode(int m);
    private:
        int mode;   // What to return when PsiSqr_intz is called
};

#endif  // WAVE_FUNCTION_H
