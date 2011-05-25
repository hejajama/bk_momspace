/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */
 

#include "config.hpp"
#include "amplitude.hpp"
#include "solver_force.hpp"
#include "solver_chebyshev.hpp"
#include "chebyshev_amplitude.hpp"
#include "tools.hpp"
#include "chebyshev.hpp"
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
    GENERATE_DATAFILE,      // Solve and print output to one huge file
    GENERATE_PLOTS,         // Print output to different files with constant rapidity
    LOGLOG_DERIVATIVE   // Calculate d ln N(k^2) / d ln(k^2)
};

enum METHOD
{
    BRUTEFORCE,             // BruteForceSolver
    CHEBYSHEV_SERIES        // ChebyshevAmplitudeSolver
};

enum CHEBYSHEV_MATRIX
{
    LOAD,
    SAVE
};

Amplitude* N;        

REAL minktsqr=DEFAULT_MINKTSQR;
REAL maxktsqr = DEFAULT_MAXKTSQR;
REAL ktsqr_mult = DEFAULT_KTSQR_MULTIPLIER;
REAL y = 0;
REAL miny = 0.0;
REAL maxy = 0.0;
REAL maxdatay = -1;
REAL delta_datay = 1;
REAL y_points = 10;
INITIAL_CONDITION ic=FTIPSAT;
MODE mode=GENERATE_DATAFILE;
bool kc=false;  // Kinematical constraint

METHOD method=CHEBYSHEV_SERIES;

// Parameters for BruteForceSolver
bool read_data=false;
string datafile="output";
int avg=0;
string file_prefix="output";

// Parameters for ChebyshevAmplitudeSolver
int chebyshev_degree=0;
std::string matrixfile="matrix.dat";
CHEBYSHEV_MATRIX cheb_matrix = LOAD;

std::stringstream infostr;

void GenerateDataFile();
void GeneratePlots();
void LogLogDerivative();

int main(int argc, char* argv[])
{
    // Print the cmdline args
    cout << "# ";
    for (int i=0; i<argc; i++)
        cout << argv[i] << " " ;
    cout << endl;

    gsl_set_error_handler(&ErrHandler);

    cout << "# " << NAME << " v. " << VERSION << " " << DATE << endl;

    if (argc>1)  if (string(argv[1])=="-help")
    {
        cout << "Usage: " << endl;
        cout << "-mode [MODE]: what to do, modes: GENERATE_DATA, GENERATE_PLOTS, LOGLOG_DERIVATIVE" << endl;
        cout << "-method [METHOD]: what method is used to solve BK, methods: BRUTEFORCE, CHEBYSHEV" << endl;
        cout << "-output [prefix]: set output file prefix, filenames are prefix_y[rapidity].dat" << endl;
        cout << "-miny, -maxy: rapidity values to solve" << endl;
        //cout << "-minktsqr, -maxktsqr: range of k_T^2 to plot, doesn't affect to limits when solving BK" << endl;
        cout << "-ic [initial condition]: set initial condition, possible ones are FTIPSAT, INVPOWER " << endl;
        cout << "-kc: apply kinematical constraint" << endl;
        cout << "-avg [avgs]: number or averagements" << endl;
        cout << "-data [datafile]: read data from datafiles from path datafile_y[rapdity].dat" << endl;
        cout << "  -maxdatay [yval]: set maximum y value for datafiles, -delta_datay [yval] difference of yvals for datafiles" << endl;
        cout << "-y [yval]: rapidity value for e.g. loglog derivative" << endl;
        cout << "-load_matrix [filename], -save_matrix [filename]: load/save coefficient matrix (CHEBYSHEV method)" << endl;
        return 0;
    }

    /*******************
     * Handle parameters
     ******************/

    for (int i=1; i<argc; i++)
    {
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
        else if (string(argv[i])=="-y")
			y=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxdatay")
            maxdatay = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-delta_datay")
            delta_datay = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-load_matrix")
        {
            cheb_matrix = LOAD;
            method = CHEBYSHEV_SERIES;  // Implicitly clear
            matrixfile = argv[i+1];
        }
        else if (string(argv[i])=="-save_matrix")
        {
            cheb_matrix = SAVE;
            method = CHEBYSHEV_SERIES;  // Implicitly clear
            matrixfile = argv[i+1];
        }
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
                mode=GENERATE_DATAFILE;
            else if (string(argv[i+1])=="LOGLOG_DERIVATIVE")
                mode=LOGLOG_DERIVATIVE;
			else if (string(argv[i+1])=="GENERATE_PLOTS")
				mode=GENERATE_PLOTS;
            else
            {
                cerr << "Unrecognized mode " << argv[i+1] << ", exiting..." << endl;
                return -1;
            }

        }
        else if (string(argv[i])=="-method")
        {
            if (string(argv[i+1])=="BRUTEFORCE")
                method = BRUTEFORCE;
            else if (string(argv[i+1])=="CHEBYSHEV")
                method = CHEBYSHEV_SERIES;
            else
                cerr << "Mode " << argv[i+1] << " is not valid!" << endl;


        }
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }


    /****************************************
     * Intialize Amplitude (=BK solver)
     ***************************************/

    if (method==BRUTEFORCE) N = new BruteForceSolver;
    else if (method==CHEBYSHEV_SERIES) N = new ChebyshevAmplitudeSolver;

    N->SetInitialCondition(ic);
    N->SetKinematicConstraint(kc);
    N->SetMaxY(maxy);
    N->SetMaxKtsqr(maxktsqr);
    if (N->YPoints()<2 and method==BRUTEFORCE and mode==GENERATE_DATAFILE)
    {
        cerr << "There must be at least 3 ypoints to evaluate " << endl;
        return -1;
    }
    N->Initialize();

    int ktsqrpoints = N->KtsqrPoints();

    infostr << "# Yrange [" << miny << ", " << maxy << "], k_T^2 limits for "
        << "output are [" << minktsqr << ", " << maxktsqr << "]" << endl;

    infostr << "# Kinematical constraint is ";
    if (kc) infostr << "applied"; else infostr << "not applied"; infostr << endl;

    infostr << "# Initial condition: " << N->InitialConditionStr() <<  endl;
    infostr << "# Grid size: ktsqrpoints x ypoints = " << N->KtsqrPoints() << " x " << N->YPoints()
        << " = " << N->KtsqrPoints()*N->YPoints() << endl;
    infostr << "# Number of averagements: " << avg << endl;
    cout << infostr.str() << endl;

    switch(method)
    {
        case BRUTEFORCE:
            N->SetNumberOfAveragements(avg);
            if (mode != GENERATE_DATAFILE)
            {
                cout << "# Reading data from file " << datafile << endl;
                ((BruteForceSolver*)N)->ReadData(datafile);    
                infostr << "# Data read from file " << datafile << endl;
            }
            else
            {
                cout << "# Generating data..." << endl;
                ((BruteForceSolver*)N)->Solve(maxy);
            }
            
            break;
        case CHEBYSHEV_SERIES:
            ((ChebyshevAmplitudeSolver*)N)->SetChebyshevDegree(200);
            ((ChebyshevAmplitudeSolver*)N)->SetBoundaryCondition(CHEBYSHEV_ZERO);
            ((ChebyshevAmplitudeSolver*)N)->Prepare();

            if (cheb_matrix==LOAD)
            {
                cout << "# Loading coefficient matrix from file " << matrixfile << endl;
                ((ChebyshevAmplitudeSolver*)N)->LoadMatrix(matrixfile);
            }
            else
            {
                cout << "# Solving coefficient matrix and saving it in file " << matrixfile << endl;
                ((ChebyshevAmplitudeSolver*)N)->SolveMatrix();
                ((ChebyshevAmplitudeSolver*)N)->SaveMatrix(matrixfile);
                cout << "# Coefficient matrix saved in file " << matrixfile << endl;
            }
            cout << "# Evolving in rapidity" << endl;
            ((ChebyshevAmplitudeSolver*)N)->Solve(maxy);
            
            
            break;
    }

    /******
     * Do some science with results
     *****/
    

    if (mode == GENERATE_DATAFILE)
        GenerateDataFile();
    else if (mode==GENERATE_PLOTS)
        GeneratePlots();
    else if (mode==LOGLOG_DERIVATIVE)
        LogLogDerivative();



    delete N;
    return 0;
}


/*
 * Generate plots, can use any kind of solver
 * (depends only on methods avaiable via Amplitude class
 */
void GeneratePlots()
{
    cout << "Saving data to files " << file_prefix << "_y.dat" << endl;
    int ktsqrpoints = N->KtsqrPoints();
    for (int yind=0; yind <= y_points; yind++)
    {
        REAL tmpy = miny + (maxy-miny)/(REAL)(y_points) * yind;
        std::ofstream output;
        std::stringstream s;
        s << file_prefix << "_y" << tmpy << ".dat";
        string fname; s >> fname;
        output.open(fname.c_str());
        output << infostr.str();

        for (int i=0; i<ktsqrpoints-1; i++)
		{
			REAL tmpktsqr = minktsqr*std::pow(ktsqr_mult, i);
			output << N->N(tmpktsqr, tmpy) << endl;
        }
        output.close();
    }
}

/*
 * Generate datafile
 * Used usually with BruteForceSolver, but works with any Amplitude-derived class
 */
void GenerateDataFile()
{
    int ktsqrpoints = N->KtsqrPoints();
    std::ofstream output;
    std::stringstream s; s << file_prefix; s << ".dat";
    output.open(s.str().c_str());
    output << infostr.str();
    cout << "Saving data to file " << s.str() << endl;
    for (int yind=0; yind <= y_points; yind++)
    {
        REAL tmpy = miny + (maxy-miny)/(REAL)(y_points) * yind;

        
        output << "###" << tmpy << endl;
        output << "# ktsqr    N(ktsqr, y=" << tmpy << ")" << endl;
		   
			
        for (int i=0; i<ktsqrpoints-1; i++)
		{
			REAL tmpktsqr = minktsqr*std::pow(ktsqr_mult, i);
			output << N->N(tmpktsqr, tmpy) << endl;
        }
    }
    output.close();
}

/*
 * LogLog-derivative
 * Calculates d ln(N(k^2)) / d ln(k^2) at given y
 * Doesn't depend on the solver used
 */
void LogLogDerivative()
{
    int ktsqrpoints = N->KtsqrPoints();
    cout << "# d ln N(k^2) / d ln (k^2), y=" << y << endl;
	cout << "#ktsqr derivative " << endl;
		
	for (int i=0; i<ktsqrpoints-1; i++)
	{
		REAL tmpktsqr = minktsqr*std::pow(ktsqr_mult, i);
		cout << tmpktsqr << " " << N->LogLogDerivative(tmpktsqr, y) << endl;
	}
}
