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

        // New rapidity?
        if (line.substr(0,3)=="###")
        {
            if (tmpvec.size()>0)
                data.push_back(tmpvec);

            if (tmpvec.size()>0 and tmpvec.size() != ktsqrpoints)
                {
                    cerr << "File " << fname << ": read " << tmpvec.size() << " ktsqrpoints, but "
                    << "there should have been " << ktsqrpoints << " points, y=" 
                    << ". " << LINEINFO << endl;
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

void DataFile::GetData(std::vector< std::vector<REAL> > &n)
{
	/*if (n[0].size()>1)	// Contains more that intial condition (which will be overrided also)
		cerr << "Got non-empty table for amplitude values? " << LINEINFO << endl;
	*/
	n.clear();
    // Return vector where indexes are vec[ktsqr][y]
    for (unsigned int k=0; k<ktsqrpoints; k++)
    {
        std::vector<REAL> tmpvec;
        for (unsigned int y=0; y<data.size(); y++)
        {
            tmpvec.push_back(data[y][k]);
        }
        n.push_back(tmpvec);
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

REAL DataFile::DeltaY()
{
	//NOTE: Assumes that y[n+1]-y[n]=const
	return yvals[1]-yvals[0];
}
