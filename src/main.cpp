/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */
 

#include "config.hpp"
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

enum MODE
{
    GENERATE_DATA,      // Solve and print output to files
    LOGLOG_DERIVATIVE   // Calculate d ln N(k^2) / d ln(k^2)
};

int main(int argc, char* argv[])
{
    Amplitude N;
    REAL minktsqr=DEFAULT_MINKTSQR;
    REAL maxktsqr = DEFAULT_MAXKTSQR;
    REAL ktsqr_mult = DEFAULT_KTSQR_MULTIPLIER;
    REAL delta_y = DEFAULT_DELTA_Y;
    REAL y = 0;
    REAL miny = 0.0;
    REAL maxy = 1.0;
    REAL maxdatay = -1;
    REAL delta_datay = 1;
    REAL y_points = 10;
    INITIAL_CONDITION ic=FTIPSAT;
    MODE mode=GENERATE_DATA;
    bool kc=false;  // Kinematical constraint
    bool read_data=false;
    string datafile;
    int avg=0;

    gsl_set_error_handler(&ErrHandler);

    OUTPUT output = OUTPUT_FILE;
    string file_prefix="output";

    cout << "# " << NAME << " v. " << VERSION << " " << DATE << endl;

    if (argc>1)  if (string(argv[1])=="-help")
    {
        cout << "Usage: " << endl;
        cout << "-mode [MODE]: what to do, modes: GENERATE_DATA, LOGLOG_DERIVATIVE" << endl;
        cout << "-output [prefix]: set output file prefix, filenames are prefix_y[rapidity].dat" << endl;
        cout << "-miny, -maxy: rapidity values to solve" << endl;
        //cout << "-minktsqr, -maxktsqr: range of k_T^2 to plot, doesn't affect to limits when solving BK" << endl;
        cout << "-ic [initial condition]: set initial condition, possible ones are FTIPSAT, INVPOWER " << endl;
        cout << "-kc: apply kinematical constraint" << endl;
        cout << "-avg [avgs]: number or averagements" << endl;
        cout << "-data [datafile]: read data from datafiles from path datafile_y[rapdity].dat" << endl;
        cout << "  -maxdatay [yval]: set maximum y value for datafiles, -delta_datay [yval] difference of yvals for datafiles" << endl;
        return 0;
    } 

    for (int i=1; i<argc; i++)
    {
        /*if (string(argv[i])=="-minktsqr")
            minktsqr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxktsqr")
            maxktsqr = StrToReal(argv[i+1]);
        else  if (string(argv[i])=="-y")
            y = StrToReal(argv[i+1]);
        else */
        if (string(argv[i])=="-miny")
            miny = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxy")
            maxy = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-N_y")
            y_points = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-output")
            file_prefix=argv[i+1];
        else if (string(argv[i])=="-kc")
            kc=true;
        else if (string(argv[i])=="-avg")
            avg=StrToInt(argv[i+1]);
        else if (string(argv[i])=="-data")
        {
            read_data=true;
            datafile= argv[i+1];
        }
        else if (string(argv[i])=="-maxdatay")
            maxdatay = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-delta_datay")
            delta_datay = StrToReal(argv[i+1]);
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
        else if (string(argv[i])=="-mode")
        {
            if (string(argv[i+1])=="GENERATE_DATA")
                mode=GENERATE_DATA;
            else if (string(argv[i+1])=="LOGLOG_DERIVATIVE")
                mode=LOGLOG_DERIVATIVE;
            else
            {
                cerr << "Unrecognized mode " << argv[i+1] << ", exiting..." << endl;
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
    N.SetNumberOfAveragements(avg);
    N.SetMaxY(maxy);
    if (N.YPoints()<2)
    {
        cerr << "There must be at least 3 ypoints to evaluate " << endl;
        return -1;
    }
    N.Initialize();

    std::stringstream infostr;
    infostr << "# Yrange [" << miny << ", " << maxy << "], k_T^2 limits for "
        << "output are [" << minktsqr << ", " << maxktsqr << "]" << endl;
    infostr << "# Kinematical constraint is ";
    if (kc) infostr << "applied"; else infostr << "not applied"; infostr << endl;
    infostr << "# Initial condition: " << N.InitialConditionStr() <<  endl;
    infostr << "# Grid size: ktsqrpoints x ypoints = " << N.KtsqrPoints() << " x " << N.YPoints()
        << " = " << N.KtsqrPoints()*N.YPoints() << endl;
    infostr << "# Number of averagements: " << avg << endl;
    cout << infostr.str();

    
    if (mode==GENERATE_DATA)
    {
        // First 3 non-comment lines, see README for syntax reference
        infostr << "###" <<minktsqr << endl << "###" << N.KtsqrMultiplier() << endl
            << "###" << N.KtsqrPoints()-1 << endl;
        N.Solve(maxy);
    }
    else
    {
        // Read data from files


    }
    int ktsqrpoints = N.KtsqrPoints();


    for (int yind=0; yind <= y_points; yind++)
    {
        REAL tmpy = miny + (maxy-miny)/(REAL)(y_points) * yind;

        std::ofstream out;
        if (output == OUTPUT_FILE)
        {
            std::stringstream s;
            s << file_prefix << "_y" << tmpy << ".dat";
            string fname; s >> fname;
            out.open(fname.c_str());

            
            
            out << infostr.str();
            out << "# ktsqr    N(ktsqr, y=" << tmpy << ")" << endl;
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
