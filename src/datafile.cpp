/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "tools.hpp"
#include "config.hpp"
#include "datafile.hpp"
#include <fstream>
#include <sstream>
using std::ifstream;
using std::getline;
using std::stringstream;

DataFile::DataFile(string fname)
{
    filename=fname;
    ifstream file(fname.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << fname << endl;
        return;
    }
    int confid=0;
    while(!file.eof() and confid < 3)
    {
        string line;
        getline(file, line);
        if (line.substr(0, 3)=="###")
        {                    
            switch (confid)
            {
                case 0:
                    minktsqr = StrToReal(line.substr(3,line.length()-3));
                    break;
                case 1:
                    ktsqr_multiplier = StrToReal(line.substr(3,line.length()-3));
                    break;
                case 2:
                    ktsqrpoints = StrToInt(line.substr(3,line.length()-3));
                    break;
                default:
                    cerr << "File " << fname << " is formatted incorrectly!" << endl;
                    break;
            }
            confid++;
        }
    }

    // Ok, configurations are read, then read all yvals
    REAL y=-1;
    std::vector<REAL> tmpvec;
    while (!file.eof())
    {
        string line;
        getline(file, line);
        if (line.substr(0,1)=="#")
            continue;   // Comment

        // New rapidity?
        if (line.substr(0,3)=="###")
        {
            if (tmpvec.size()>0)
                data.push_back(tmpvec);

            if (tmpvec.size()>0 and tmpvec.size() != ktsqrpoints)
                {
                    cerr << "File " << fname << ": read " << tmpvec.size() << " ktsqrpoints, but "
                    << "there should have been " << ktsqrpoints << " points, y=" << y << endl;
                }

            y = StrToReal(line.substr(3,line.length()-3));
            yvals.push_back(y);
            tmpvec.clear();
            continue;   // Next line is probably new amplitude value
        }

        // Ok, so this a new amplitude value
        tmpvec.push_back(StrToReal(line));
    }

    if (data.size() != ktsqrpoints)
    {
        cerr << "File " << fname << ": read " << data.size() << " ktsqrpoints, but "
        << "there should have been " << ktsqrpoints << " points???" << endl;
    }
}

std::vector< std::vector<REAL> >& DataFile::GetData()
{
    // Return vector where indexes are vec[ktsqr][y]
    std::vector< std::vector<REAL> > *result = new std::vector< std::vector<REAL> >;
    for (int k=0; k<ktsqrpoints; k++)
    {
        std::vector<REAL> tmpvec;
        for (int y=0; y<data.size(); y++)
        {
            tmpvec.push_back(data[y][k]);
        }
        (*result).push_back(tmpvec);
    }
    return (*result);
}

REAL DataFile::MinKtsqr()
{
    return minktsqr;
}

REAL DataFile::KtsqrMultiplier()
{
    return ktsqr_multiplier;
}

REAL DataFile::KtsqrPoints()
{
        return ktsqrpoints;
}
