/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _DATAFILE_H
#define _DATAFILE_H

#include <tools/tools.hpp>
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
		REAL MaxY();

        void GetData(std::vector< std::vector<REAL> > &ln_n,
            std::vector<REAL> &rapidities);

    private:
        string filename;
        std::vector<std::vector <REAL> > data;
        std::vector<REAL> yvals;
        REAL minktsqr;
        REAL ktsqr_multiplier;
        unsigned int ktsqrpoints;

};

#endif
