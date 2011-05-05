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
    while(!file.eof())
    {
        string line;
        getline(file, line);
        if (line.substr(0, 3)=="###")
        {                    
            switch (confid)
            {
                case 0:
                    minktsqr = StrToReal(line);
                    break;
                case 1:
                    ktsqr_multiplier = StrToReal(line);
                    break;
                case 2:
                    ktsqrpoints = StrToInt(line);
                    break;
                default:
                    cerr << "File " << fname << " is formatted incorrectly!" << endl;
                    break;
            }
            confid++;
        }

        if (line.substr(0,1)=="#")
            continue;   // Comment

        // Ok, so this is ktsqr N pair
        stringstream ss; string tmp;
        ss << line;
        ss >> tmp; tmp=""; ss >> tmp;   // Second word is what we are interested in
        data.push_back(StrToReal(tmp));
    }

    if (data.size() != ktsqrpoints)
    {
        cerr << "File " << fname << ": read " << data.size() << " ktsqrpoints, but "
        << "there should have been " << ktsqrpoints << " points???" << endl;
    }
}

std::vector<REAL>& DataFile::GetData()
{
    return data;
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
