/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _AMPLITUDE_H
#define _AMPLITUDE_H

#include "config.hpp"
#include <vector>
#include <cmath>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
//#include <bci.h>

const REAL DEFAULT_DELTA_Y=0.1;
const unsigned int DEFAULT_MAXY=50;
const REAL DEFAULT_KTSQR_MULTIPLIER = 1.01;  // ktsqr_{i+1} = ktsqr_i*KTSQR_MULTIPLIER
//const REAL DEFAULT_KTSQR_MULTIPLIER = 1.05;
const REAL DEFAULT_MINKTSQR = 1e-8; // orig: 1e-8
const REAL DEFAULT_MAXKTSQR = 1e10; // orig: 1e10

//const unsigned int POINTS_Y= (int)(MAXY/DELTA_Y);
//const unsigned int POINTS_KTSQR = (int)(std::log(MAXKTSQR/MINKTSQR) / std::log(KTSQR_MULTIPLIER) );

const REAL KTSQRINTACCURACY = 0.005;
const int KTSQRINTITERATIONS = 9000;
const int INTERPOLATION_POINTS = 10;
const int INTERPOLATION_POINTS_DER=200;

enum INITIAL_CONDITION
{
    FTIPSAT,    // \int d^2r 1/(2\pi r^2) exp(ik.r) 2(1-exp(-r^2))
    INVPOWER,   // 1/(k^2 + 1), as in BK in full mom. space, hep-ph/0504080
    INVPOWER4,  // 1/(k^4+1), very arbitrary
    GAUSS       // Exp[-(Log[k^2] + 2)^2/5], hep-ph/0110325
};
    

class Amplitude
{
    public:
        Amplitude();
        ~Amplitude();
        void Initialize();
		void Clear();

        // Notice: N is virtual, so subclasses may use their own methods
        // to compute the amplitude (e.g. ChebyshevAmplitudeSolver), and
        // in that case they may not use n[][]-table or AddDataPoint-routines
        // at all!
        virtual REAL N(REAL ktsqr, REAL y, bool bspline=false);
        void IntializeBSpline(int ktsqrind, REAL rapidity);
        REAL BSplineAmplitude(REAL ktsqr, REAL rapidity);
        //REAL BSplineAmplitude(REAL ktsqr, REAL *ktsqrarray, REAL *narray, uint points);
        void AddDataPoint(int ktsqrindex, int yindex, REAL value, REAL der);
        
        REAL LogLogDerivative(REAL ktsqr, REAL y);

        int ReadData(string file);
        

        virtual void Solve(REAL maxy)=0;

        REAL Ktsqrval(unsigned int i);
        REAL Yval(unsigned int i);
        void SetInitialCondition(INITIAL_CONDITION i);
        void SetKinematicConstraint(bool kc);
        bool KinematicalConstraint();
        string InitialConditionStr();
        void SetNumberOfAveragements(int avg);

        void SetMinKtsqr(REAL mkt);
        void SetMaxKtsqr(REAL mkt);
        void SetKtsqrMultiplier(REAL m);
        void SetMaxY(REAL y);
        void SetDeltaY(REAL dy);

        void SetRunningCoupling(bool rc);
        bool RunningCoupling();

        unsigned int YPoints();
        unsigned int KtsqrPoints();
        REAL KtsqrMultiplier();
        REAL DeltaY();
        virtual REAL MinKtsqr();
        virtual REAL MaxKtsqr();

        REAL InitialCondition(REAL ktsqr);  // N() at y=0
        
    protected:

        REAL BSplineDerivative(REAL ktsqr, REAL* ktsqrarray, REAL* narray, uint points);
    
        // n[ktsqr][y]
        std::vector< std::vector<REAL> > n;

        std::vector<REAL> ktsqrvals;
        std::vector<REAL> yvals;

        // derivatives
        std::vector< std::vector<REAL> > derivatives;

        INITIAL_CONDITION ic;
        bool kinematic_constraint;
        int averages;

        REAL minktsqr;
        REAL maxktsqr;
        REAL ktsqr_multiplier;
        REAL maxy;
        REAL delta_y;

        bool running_coupling;

        bool datafile;	// True if data is read from external file

        bool adams_method;   // Whether the Adams method should be used when solvin DE

        // Bspline interpolation
        REAL bspline_y;
        gsl_bspline_workspace *bw;
        gsl_vector *B;
        gsl_vector *c;
        gsl_matrix *X;
        gsl_matrix *cov;
        gsl_multifit_linear_workspace *mw;
        int interpolation_start;    // These values are used to cache interpolation
        int interpolation_end;
        REAL interpolation_rapidity;
        
};

#endif
