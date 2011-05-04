/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */
 

#include "config.hpp"
#include "function.hpp"
#include "amplitude.hpp"
#include "tools.hpp"
#include <gsl/gsl_errno.h>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>

using std::string;

const string VERSION = "0.1-dev";
const string DATE = "2011-xx-xx";
const string NAME = "BK";

enum OUTPUT
{
    OUTPUT_FILE, OUTPUT_TERMINAL
};

int main(int argc, char* argv[])
{
    Amplitude N;
    REAL minktsqr=MINKTSQR;
    REAL maxktsqr = MAXKTSQR;
    REAL ktsqr_mult = KTSQR_MULTIPLIER;
    REAL y = 0;
    REAL miny = 0.0;
    REAL maxy = 1.0;
    REAL y_points = 10;
    INITIAL_CONDITION ic=FTIPSAT;
    bool kc=false;  // Kinematical constraint

    gsl_set_error_handler(&ErrHandler);

    OUTPUT output = OUTPUT_FILE;
    string file_prefix="output_";

    cout << "# " << NAME << " v. " << VERSION << " " << DATE << endl;

    if (argc>1)  if (string(argv[1])=="-help")
    {
        cout << "Usage: " << endl;
        cout << "-output [prefix]: set output file prefix, filenames are prefix_y[rapidity].dat" << endl;
        cout << "-miny, -maxy: rapidity values to solve" << endl;
        cout << "-minktsqr, -maxktsqr: range of k_T^2 to plot, doesn't affect to limits when solving BK" << endl;
        cout << "-ic [initial condition]: set initial condition, possible ones are FTIPSAT, INVPOWER " << endl;
        cout << "-kc: apply kinematical constraint" << endl;
        return 0;
    } 

    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-minktsqr")
            minktsqr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxktsqr")
            maxktsqr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-y")
            y = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-miny")
            miny = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxy")
            maxy = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-N_y")
            y_points = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-output")
            file_prefix=argv[i+1];
        else if (string(argv[i])=="-kc")
            kc=true;
        else if (string(argv[i])=="-ic")
        {
            if (string(argv[i+1])=="FTIPSAT")
                ic = FTIPSAT;
            else if (string(argv[i+1])=="INVPOWER")
                ic = INVPOWER;
            else
            {
                cerr << "Unrecognized initial condition " << argv[i+1] << ", exiting..." << endl;
                return -1;
            }
        }
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }
    N.SetInitialCondition(ic);
    N.SetKinematicConstraint(kc);
    N.Initialize();

    std::stringstream infostr;
    infostr << "# Solving BK using yrange [" << miny << ", " << maxy << "], k_T^2 limits for "
        << "output are [" << minktsqr << ", " << maxktsqr << "]" << endl;
    infostr << "# Kinematical constraint is ";
    if (kc) infostr << "applied"; else infostr << "not applied"; infostr << endl;
    infostr << "# Initial condition: " << N.InitialConditionStr() <<  endl;
    infostr << "# Grid size: ktsqrpoints x ypoints = " << POINTS_KTSQR << " x " << (int)(maxy/DELTA_Y)
        << " = " << POINTS_KTSQR*(int)(maxy/DELTA_Y) << endl;
    cout << infostr.str();

    //cout << N.Ktsqrval(1000);
   //N.RapidityDerivative(20.0, 0); return 0;
    /*for (REAL ktsqr=MINKTSQR; ktsqr<=MAXKTSQR; ktsqr*=KTSQR_MULTIPLIER )
    {
        cout << ktsqr << " " << N.RapidityDerivative(ktsqr, 0.0) << endl;
    }return 0;*/


    N.Solve(maxy);
/*
    for (REAL ktsqr=MINKTSQR; ktsqr<=MAXKTSQR; ktsqr*=KTSQR_MULTIPLIER )
    {
        cout << ktsqr << " " << N.RapidityDerivative(ktsqr, 0.1) << endl;
    }return 0;
  */
//return 0;
    int ktsqrpoints = static_cast<int>( std::log(maxktsqr/minktsqr) / std::log(ktsqr_mult) );


    for (int yind=0; yind <= y_points; yind++)
    {
        REAL tmpy = miny + (maxy-miny)/(REAL)(y_points) * yind;

        std::ofstream out;
        if (output == OUTPUT_FILE)
        {
            std::stringstream s;
            s << file_prefix << "y" << tmpy << ".dat";
            string fname; s >> fname;
            out.open(fname.c_str());
            out << infostr.str();
            out << "# ktsqr    N(ktsqr, y=" << tmpy << endl;
        }
        
       // cout << endl << endl << endl;
        
        for (int i=0; i<ktsqrpoints-1; i++)
        {
            REAL tmpktsqr = minktsqr*std::pow(ktsqr_mult, i);
            if (output == OUTPUT_TERMINAL)
            {
                cout << tmpktsqr << " " << N.N(tmpktsqr, tmpy) << endl;
            } else if (output == OUTPUT_FILE)
            {
                out << tmpktsqr << " " << N.N(tmpktsqr, tmpy) << endl;
            }
        }
        out.close();
    }
/*
    N.AddDataPoint(0, 1, 8);
    N.AddDataPoint(0, 2, 3);
    N.AddDataPoint(1, 1, 6);
    //cout << N.N(0,0) << " "
    //    << N.N(0.05, DELTA_Y*0.5) << endl;
    cout << endl << N.RapidityDerivative(1,0) << endl;
  */  

    return 0;
}
