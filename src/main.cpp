/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */
 

#include "config.hpp"
#include "amplitude.hpp"
#include "solver_force.hpp"
#include "solver_force2.hpp"
#include "solver_chebyshev.hpp"
#include "chebyshev_amplitude.hpp"
#include "tools.hpp"
#include "chebyshev.hpp"
#include "hankel.hpp"
#include "interpolation.hpp"
#include <gsl/gsl_errno.h>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
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
    GENERATE_SINGLE_PLOT,   // Print amplitude at a given rapidity
    GENERATE_SINGLE_RPLOT,  // Print amplitude in r-space at a given rapidity
    LOGLOG_DERIVATIVE,   // Calculate d ln N(k^2) / d ln(k^2)
    SATURATION_SCALE        // Calculate Q_s as a function of y up to maxy
};

enum METHOD
{
    BRUTEFORCE,             // BruteForceSolver
    BRUTEFORCE2,            // BruteForceSolver2
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
uint ktsqrpoints = static_cast<uint>(std::log(maxktsqr/minktsqr) / std::log(ktsqr_mult));
REAL y = 0;
REAL miny = 0.0;
REAL maxy = 0.0;
REAL maxdatay = -1;
REAL delta_datay = 1;
REAL y_points = 10;
INITIAL_CONDITION ic=FTIPSAT;
MODE mode=GENERATE_DATAFILE;
bool kc=false;  // Kinematical constraint
RUNNING_COUPLING running_coupling=CONSTANT;
bool scale_sat=false;  // Scale k_T with the scaturation scale

METHOD method=CHEBYSHEV_SERIES;

// Parameters for BruteForceSolver
bool read_data=false;
string datafile="output";
int avg=0;
string file_prefix="output";
bool rungekutta = false;
bool adams = false;

// Parameters for ChebyshevAmplitudeSolver
int chebyshev_degree=100;
std::string matrixfile="matrix.dat";
CHEBYSHEV_MATRIX cheb_matrix = LOAD;

std::stringstream infostr;

void GenerateDataFile();
void GeneratePlots();
void LogLogDerivative();
void SinglePlot();
void SinglePlotR();     // FT to r-space
void Clear();
void SaturationScale();

int main(int argc, char* argv[])
{    
    // Print the cmdline args
    std::stringstream cmdline;
    cmdline << "# ";
    for (int i=0; i<argc; i++)
        cmdline << argv[i] << " " ;
    cmdline << endl;
    cout << cmdline.str();

    gsl_set_error_handler(&ErrHandler);

    cout << "# " << NAME << " v. " << VERSION << " " << DATE << endl;

    if (argc>1)  if (string(argv[1])=="-help")
    {
        cout << "Usage: " << endl;
        cout << "-mode [MODE]: what to do, modes: GENERATE_DATA, GENERATE_PLOTS, SINGLE_PLOT, SINGLE_RPLOT, LOGLOG_DERIVATIVE " << endl;
        cout << "              SATURATION_SCALE" << endl;
        cout << "-method [METHOD]: what method is used to solve BK, methods: BRUTEFORCE[2], CHEBYSHEV" << endl;
        cout << "-output [prefix]: set output file prefix, filenames are prefix_y[rapidity].dat" << endl;
        cout << "-miny, -maxy: rapidity values to solve" << endl;
        //cout << "-minktsqr, -maxktsqr: range of k_T^2 to plot, doesn't affect to limits when solving BK" << endl;
        cout << "-ic [initial condition]: set initial condition, possible ones are FTIPSAT, INVPOWER, INVPOWER4, GAUSS " << endl;
        cout << "-kc: apply kinematical constraint" << endl;
        cout << "-rc [METHOD]: apply running coupling, methods: CONSTANT, PARENT, MINK, MAXK" << endl;
        cout << "-adams_method: use Adams method with BRUTEFORCE)" << endl;
        cout << "-avg [avgs]: number or averagements" << endl;
        cout << "-data [datafile]: read data from datafiles from path datafile_y[rapdity].dat" << endl;
        cout << "  -maxdatay [yval]: set maximum y value for datafiles, -delta_datay [yval] difference of yvals for datafiles" << endl;
        cout << "-y [yval]: rapidity value for e.g. loglog derivative" << endl;
        cout << "-load_matrix [filename], -save_matrix [filename]: load/save coefficient matrix (CHEBYSHEV method)" << endl;
        cout << "-chebyshev_degree [number]: number of basis vectors" << endl;
        cout << "-minktsqr [val], -maxktsqr [val], -ktsqrpoints [val]" << endl;
        cout << "-scale_sat: scale k_T by saturation scale" << endl;
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
        else if (string(argv[i])=="-rc")
        {
            if (string(argv[i+1])=="CONSTANT")
                running_coupling = CONSTANT;
            else if (string(argv[i+1])=="PARENT")
                running_coupling = PARENT_DIPOLE;
            else if (string(argv[i+1])=="MAXK")
                running_coupling = MAXK;
            else if (string(argv[i+1])=="MINK")
                running_coupling = MINK;
            else
            {
                cerr << "Unrecognized running coupling " << argv[i+1] << ", exiting..." << endl;
                return -1;
            }
        }
        else if (string(argv[i])=="-adams_method")
            adams = true;
        else if (string(argv[i])=="-scale_sat")
            scale_sat=true;
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
        else if (string(argv[i])=="-minktsqr")
            minktsqr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxktsqr")
            maxktsqr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-ktsqrpoints")
            ktsqrpoints = StrToInt(argv[i+1]);
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
        else if (string(argv[i])=="-chebyshev_degree")
            chebyshev_degree = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-ic")
        {
            if (string(argv[i+1])=="FTIPSAT")
                ic = FTIPSAT;
            else if (string(argv[i+1])=="INVPOWER")
                ic = INVPOWER;
            else if (string(argv[i+1])=="INVPOWER4")
                ic = INVPOWER4;
            else if (string(argv[i+1])=="GAUSS")
                ic = GAUSS;
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
            else if (string(argv[i+1])=="SINGLE_PLOT")
                mode=GENERATE_SINGLE_PLOT;
            else if (string(argv[i+1])=="SINGLE_RPLOT")
                mode=GENERATE_SINGLE_RPLOT;
            else if (string(argv[i+1])=="SATURATION_SCALE")
                mode = SATURATION_SCALE;
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
            else if (string(argv[i+1])=="BRUTEFORCE2")
                method = BRUTEFORCE2;
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
    else if (method == BRUTEFORCE2) N = new BruteForceSolver2;
    else if (method==CHEBYSHEV_SERIES) N = new ChebyshevAmplitudeSolver;

    N->SetInitialCondition(ic);
    N->SetRunningCoupling(running_coupling);
    N->SetKinematicConstraint(kc);
    N->SetMaxY(maxy);
    N->SetMaxKtsqr(maxktsqr);
    N->SetMinKtsqr(minktsqr);
    N->SetKtsqrMultiplier( std::pow( maxktsqr/minktsqr,
        1.0/static_cast<REAL>(ktsqrpoints) ) );
    
    if (N->YPoints()<2 and method==BRUTEFORCE and mode==GENERATE_DATAFILE)
    {
        cerr << "There must be at least 3 ypoints to evaluate " << endl;
        return -1;
    }
    N->Initialize();

    ktsqrpoints = N->KtsqrPoints();

    infostr << cmdline.str();
    infostr << "# Yrange [" << miny << ", " << maxy << "], k_T^2 limits for "
        << "output are [" << minktsqr << ", " << maxktsqr << "]" << endl;

    infostr << "# Kinematical constraint is ";
    if (kc) infostr << "applied"; else infostr << "not applied"; infostr << endl;
    infostr << "# Running coupling is " << N->RunningCouplingStr() << endl;

    infostr << "# Initial condition: " << N->InitialConditionStr() <<  endl;
    infostr << "# Grid size: ktsqrpoints x ypoints = " << N->KtsqrPoints() << " x " << N->YPoints()
        << " = " << N->KtsqrPoints()*N->YPoints() << endl;
    infostr << "# Number of averagements: " << avg << endl;
    cout << infostr.str() << endl;

    if (y>0)    // We should evolve amplitude up to y
        maxy = y + DEFAULT_DELTA_Y; 

    switch(method)
    {
        case BRUTEFORCE2:
            cerr << "BRUTEFORCE2, shoudn't be here!" << endl;
            N->SetNumberOfAveragements(avg);
            if (mode != GENERATE_DATAFILE)
            {
                cout << "# Reading data from file " << datafile << endl;
                ((BruteForceSolver2*)N)->ReadData(datafile);
                minktsqr = N->MinKtsqr();
                maxktsqr = N->MaxKtsqr();
				ktsqr_mult = N->KtsqrMultiplier();
                ktsqrpoints = N->KtsqrPoints();
                cout << "# Data read from file " << datafile <<
                    " ktsqrlimits " << minktsqr << " - " << maxktsqr
                    << ", points " << ktsqrpoints << ", multiplier "
                    << N->KtsqrMultiplier() << endl;
            }
            else
            {
                cout << "# Generating data..." << endl;
                ((BruteForceSolver2*)N)->Solve(maxy);
            }
            
            break;
        case BRUTEFORCE:
            N->SetNumberOfAveragements(avg);
            if (mode != GENERATE_DATAFILE)
            {
                cout << "# Reading data from file " << datafile << endl;
                ((BruteForceSolver*)N)->ReadData(datafile);
                minktsqr = N->MinKtsqr();
                maxktsqr = N->MaxKtsqr();
				ktsqr_mult = N->KtsqrMultiplier();
                ktsqrpoints = N->KtsqrPoints();
                cout << "# Data read from file " << datafile <<
                    " ktsqrlimits " << minktsqr << " - " << maxktsqr
                    << ", points " << ktsqrpoints << ", multiplier "
                    << N->KtsqrMultiplier() << endl;
            }
            else
            {
                if (adams)
                {
                    ((BruteForceSolver*)N)->SetAdamsMethod(true);
                    infostr << "# Adams methdo is used" << endl;
                    cout << "# Adams methdo is used" << endl;
                }
                cout << "# Generating data..." << endl;
                ((BruteForceSolver*)N)->Solve(maxy);
            }
            
            break;
        case CHEBYSHEV_SERIES:
            ((ChebyshevAmplitudeSolver*)N)->SetBoundaryCondition(CHEBYSHEV_ZERO);

            if (cheb_matrix==LOAD)
            {
                cout << "# Loading coefficient matrix from file " << matrixfile << endl;
                ((ChebyshevAmplitudeSolver*)N)->LoadMatrix(matrixfile);
                cout << "# Matrix loaded, M=" << ((ChebyshevAmplitudeSolver*)N)->M() << endl;
                cout << "# Intializing ChebyshevAmplitudeSolver environment..." << endl;
                maxktsqr = N->MaxKtsqr();
                minktsqr = N->MinKtsqr();
                ((ChebyshevAmplitudeSolver*)N)->Prepare();

                // Check if we are asked to use fewer number of cheb polynomials
                // as it would be possible
                if (chebyshev_degree < ((ChebyshevAmplitudeSolver*)N)->ChebyshevDegree())
                    ((ChebyshevAmplitudeSolver*)N)->SetChebyshevDegree(chebyshev_degree);
            }
            else
            {
                cout << "# Intializing ChebyshevAmplitudeSolver environment..." << endl;
                ((ChebyshevAmplitudeSolver*)N)->SetChebyshevDegree(chebyshev_degree);
                ((ChebyshevAmplitudeSolver*)N)->Prepare();
                cout << "# Solving coefficient matrix and saving it in file " << matrixfile << endl;
                ((ChebyshevAmplitudeSolver*)N)->SolveMatrix();
                ((ChebyshevAmplitudeSolver*)N)->SaveMatrix(matrixfile);
                cout << "# Coefficient matrix saved in file " << matrixfile << endl;
                Clear();
                return 0;
            }
            cout << "# Chebyshev degree is " << ((ChebyshevAmplitudeSolver*)N)->ChebyshevDegree() << endl;
            cout << "# Evolving in rapidity" << endl;
            
            ((ChebyshevAmplitudeSolver*)N)->Solve(maxy);
            
            
            break;
    }

    /******
     * Do some science with results
     *****/



     ////////DEBUG
    /*ChebyshevAmplitudeSolver* ch = (ChebyshevAmplitudeSolver*)N;
    for (unsigned int i=0; i<=ch->ChebyshevDegree(); i++)
    {
        cout << ch->Coefficient(600, i) << endl;
    }
    return 0;
    */

     /////////////
    

    if (mode == GENERATE_DATAFILE)
        GenerateDataFile();
    else if (mode==GENERATE_PLOTS)
        GeneratePlots();
    else if (mode==LOGLOG_DERIVATIVE)
        LogLogDerivative();
    else if (mode==GENERATE_SINGLE_PLOT)
        SinglePlot();
    else if (mode==GENERATE_SINGLE_RPLOT)
        SinglePlotR();
    else if (mode==SATURATION_SCALE)
        SaturationScale();


    Clear();
    return 0;
}

void Clear()
{
    delete N;
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

        for (int i=0; i<N->KtsqrPoints()-1; i++)
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
    cout << "# Saving data to file " << s.str() << endl;
    output << "# Running coupling: " << N->RunningCouplingStr() << endl;

    // Metadata (see README for the syntax)
    output << "###" << std::scientific << std::setprecision(15) << N->MinKtsqr() << endl;
    output << "###" << std::scientific << std::setprecision(15) << N->KtsqrMultiplier() << endl;
    output << "###" << ktsqrpoints << endl;

    unsigned int yp = static_cast<int>(maxy/N->DeltaY())+1;
    
    
    for (unsigned int yind=0; yind <= yp; yind++)
    {
        //REAL tmpy = miny + (maxy-miny)/(REAL)(y_points) * yind;
        REAL tmpy = N->Yval(yind);

        
        output << "###" << tmpy << endl;
        
		   
			
        for (int i=0; i<ktsqrpoints; i++)
		{
			REAL tmpktsqr = minktsqr*std::pow(N->KtsqrMultiplier(), i);
			output << std::scientific << std::setprecision(15) << N->N(tmpktsqr, tmpy) << endl;
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
    cout << "# Saturation scale: k_T:" << endl;
    cout << "### " << N->SaturationScale(y) << endl;
	cout << "#ktsqr derivative " << endl;

	for (int i=0; i<ktsqrpoints-1; i++)
	{
		REAL tmpktsqr = minktsqr*std::pow(ktsqr_mult, i);
		cout << tmpktsqr << " " << N->LogLogDerivative(tmpktsqr, y) << endl;
	}
}

/*
 * Plot amplitude at given rapidity
 */
void SinglePlot()
{
    cout << "# y=" << y << ", ic=" << N->InitialConditionStr() << endl;
    cout << "# Saturation scale: k_T:" << endl;
    cout << "### " << N->SaturationScale(y) << endl;
    cout << "# ktsqr amplitude initial_condition bspline_amplitude" << endl;
    ktsqr_mult = N->KtsqrMultiplier();

    for (int i=0; i<N->KtsqrPoints()-1; i+=10)
    {
        REAL tmpktsqr = N->MinKtsqr()*std::pow(ktsqr_mult, i);
        cout << std::scientific << std::setprecision(15) << tmpktsqr << " " <<
            N->N(tmpktsqr, y) << " " << N->InitialCondition(tmpktsqr)
            << " " << N->N(tmpktsqr, y, true) << endl;
    }
}

/*
 * Plot amplitude in r-space at given rapidity
 */

REAL Inthelperf_ft(REAL ktqr, void* p)
{
    return 0;
}
void SinglePlotR()
{
    cout << "# y=" << y << ", ic=" << N->InitialConditionStr() << endl;
    Hankel transformed(N, y, 100);
    transformed.PrintRAmplitude();

}

/*
 * Plot Q_s as a function of y up to maxy
 */
void SaturationScale()
{
    cout << "# Saturation scale Q_s as a function of y" << endl;
    cout << "# y  Bspline-Q_s Q_s " << endl;

    // Bspline interpolation
    int points = static_cast<int>(maxy/0.1);
    REAL *yarray = new REAL[points];
    REAL *qsarray = new REAL[points];

    const int interpolation_points=20;
    int interpolation_start, interpolation_end;
    
    for (int i=0; i<points; i++)
    {
        if (i-interpolation_points/2<0)
        {
            interpolation_start=0; interpolation_end=interpolation_points-1;
        }
        else if (i+interpolation_points/2 > points)
        {
            interpolation_end = points-1;
            interpolation_start = points-1-interpolation_points;
        }
        else
        {
            interpolation_start = i-interpolation_points/2;
            interpolation_end = i+interpolation_points/2;
        }
        int interpo_points = interpolation_end - interpolation_start + 1;

        REAL *yarray = new REAL[interpo_points];
        REAL *qsarray = new REAL[interpo_points];
        for (int j=interpolation_start; j<=interpolation_end; j++)
        {
            yarray[j-interpolation_start] = j*0.1;
            qsarray[j-interpolation_start] = N->SaturationScale(j*0.1);
        }
        
        Interpolator inter(yarray, qsarray, interpo_points);
        inter.SetMethod(INTERPOLATE_BSPLINE);
        inter.Initialize();
        
        REAL tmpy = 0.1*i;
        cout << std::scientific << std::setprecision(15) << tmpy << " "
            << inter.Evaluate(tmpy) << " " << N->SaturationScale(tmpy) << endl;

        delete[] yarray;
        delete[] qsarray;
    }
}

