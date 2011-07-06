/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _SPECTRUM_H
#define _SPECTRUM_H

#include "amplitude.hpp"

/*
 * This class calculates various spectra, e.g. charge hadron multiplicity
 */

class Spectrum
{
    public:
        Spectrum(Amplitude* amp);
        // dN_ch / dydp_T^2, integrated over b and angular dep. of p_T
        REAL dNch_dydpsqr(REAL sqrts, REAL y, REAL psqr);
        // Previous integrated over p
        REAL dNch_dy(REAL sqrts, REAL y);
        
    private:
        Amplitude* N;

};



#endif
