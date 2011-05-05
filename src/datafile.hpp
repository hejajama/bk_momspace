/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _DATAFILE_H
#define _DATAFILE_H

#include "tools.hpp"
#include "config.hpp"
#include <sstream>
#include <vector>
/* Read data from datafiles, file format is defined in file README
 * NOTICE: This doesn't use the k^2 value given in each row, instead the
 * value of k^2 is computed using the fact that ktsqr(i) = ktsqr(0)*pow(multiplier, i)
 */

class DataFile
{
    public:
        DataFile(string fname);
        REAL MinKtsqr();
        REAL KtsqrMultiplier();
        REAL KtsqrPoints();

        std::vector<REAL>& GetData();

    private:
        string filename;
        std::vector<REAL> data;
        REAL minktsqr;
        REAL ktsqr_multiplier;
        unsigned int ktsqrpoints;

};

#endif
