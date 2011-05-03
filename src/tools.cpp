/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "tools.hpp"
#include "config.hpp"
#include <string>
#include <sstream>
#include <cmath>

/*
 * Str to REAL/int
 */
REAL StrToReal(std::string str)
{
    std::stringstream buff(str);
    REAL tmp;
    buff >> tmp;
    return tmp;
}

int StrToInt(std::string str)
{
    std::stringstream buff(str);
    int tmp;
    buff >> tmp;
    return tmp;
}

// GSL Error handler
int errors;
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{
    errors++;
    if (gsl_errno !=14  )  // 14 = failed to reach tolerance, handle when gsl_int 
                    // is called
        std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

/* 
 * Q^2 dependent strong coupling constant
 * Takes into account only u,d ands s quarks
 */
REAL Alpha_s(REAL Qsqr)
{
    return 12.0*M_PI/( (33.0-2.0*Nf)*log(Qsqr/LAMBDAQCD2) );
}


