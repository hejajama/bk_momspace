/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "tools.hpp"
#include "config.hpp"
#include "datafile.hpp"
#include "amplitude.hpp"
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

    //TODO: It's impossible that this condition doesn't hold
    if (confid < 3)
        cerr << "File " << fname << " doesn't have enough metadata!" << endl;

    // Ok, configurations are read, then read all yvals
    REAL y=-1;
    std::vector<REAL> tmpvec;
    while (!file.eof())
    {
        string line;
        getline(file, line);

        // New rapidity?
        if (line.substr(0,3)=="###")
        {
            if (tmpvec.size()>0)
                data.push_back(tmpvec);

            if (tmpvec.size()>0 and tmpvec.size() != ktsqrpoints)
                {
                    cerr << "File " << fname << ": read " << tmpvec.size() << " ktsqrpoints, but "
                    << "there should have been " << ktsqrpoints << " points, y=" 
                    << y << ". " << LINEINFO << endl;
                }

            y = StrToReal(line.substr(3,line.length()-3));
            yvals.push_back(y);
            tmpvec.clear();
            continue;   // Next line is probably new amplitude value
        }
        else if (line.substr(0,1)=="#")
            continue;   // Comment

        // Ok, so this a new amplitude value
        tmpvec.push_back(StrToReal(line));
    }

    if (data[0].size() != ktsqrpoints)
    {
        cerr << "File " << fname << ": read " << data.size() << " ktsqrpoints, but "
        << "there should have been " << ktsqrpoints << " points???" << LINEINFO << endl;
    }
}

void DataFile::GetData(std::vector< std::vector<REAL> > &ln_n,
                        std::vector<REAL> &rapidities)
{
	ln_n.clear();
    rapidities.clear();
    // Return vector where indexes are vec[y][ktsqr] containing ln of amplitude

    for (uint yind=0; yind < data.size(); yind++)
    {
        std::vector<REAL> tmpvec;
        for (uint kind=0; kind<ktsqrpoints; kind++)
        {
            REAL tmpln_n = std::log(data[yind][kind]);
            if (tmpln_n > MINLN_N) tmpvec.push_back(tmpln_n);
            else
                tmpvec.push_back(MINLN_N);
        }
        ln_n.push_back(tmpvec);

        rapidities.push_back(yvals[yind]);
        
    }
    
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

REAL DataFile::MaxY()
{
	return yvals[yvals.size()-1];
}



