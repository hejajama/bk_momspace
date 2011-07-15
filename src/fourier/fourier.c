/*
 * (C) notice added by Heikki Mäntysaari <heikki.mantysaari@jyu.fi>:
 * This code is downloaded from
 * http://ipht.cea.fr/Pisp/francois.gelis/Soft/Fourier/index.php
 * Although the license is not specified by the author, it is GNU GPL as
 * GSL is licensed under the GPL license which forces also this program to
 * be redistributed under the same license.
 *
 * Original author&copyright holder: Francois Gelis
 * Changes by H.M: J0zero and J1zero are now defined in separeted .h
 * files and the number of zeroes is increased
 * gsl_sum_levin_u_workspace *workspace is now intialized and destroyed
 * every time fourier_j_i is called in order to make this thread-safe
 * (but adds a bit overhead, sorry).
 */

// Last modified on June 06th, 2006 -- F. Gelis


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>

#ifdef X86
#include <fpu_control.h>
#endif

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sum.h>



typedef struct {
  double x;
#ifndef FORTRAN
  double (*function)(double,void *);
#else
  double (*function)(double *,void *);
#endif
  void *params;
} fourier_int_data;


static int Nzeros;
//static gsl_sum_levin_u_workspace *workspace;
static double epsilon1=1.0e-12;
static double epsilon=1.0e-12;


/* Table with the first 1000 zeros of J0(x) */

//double J0zero[2001]={0.0,\ ...
#include "besselj0zero.h"



//double J1zero[1001]={0.0,\ ...
#include "besselj1zero.h"


// Generic integration routine in an interval

double integrate(double a,double b,gsl_function F,\
		 gsl_integration_workspace *gsl_wksp){
  double error,term_n;

  gsl_integration_qag(&F,a,b,0.0,epsilon1,4000,6,gsl_wksp,&term_n,&error);

  return term_n;
}


// Integrands

double fourier_J0_int(double k,void *param){
  fourier_int_data *dt=(fourier_int_data*)param;

#ifndef FORTRAN
  return j0(k*(dt->x))*(dt->function)(k,dt->params);
#else
  return j0(k*(dt->x))*(dt->function)(&k,dt->params);
#endif
}

double fourier_J1_int(double k,void *param){
  fourier_int_data *dt=(fourier_int_data*)param;

#ifndef FORTRAN
  return j1(k*(dt->x))*(dt->function)(k,dt->params);
#else
  return j1(k*(dt->x))*(dt->function)(&k,dt->params);
#endif
}


// 2-d Fourier integral of the function "func", starting at the Ni-th
// zero of the Bessel function J0

double fourier_j_i(x,Ni,func,param,n)
     double x;
     int Ni;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
     int n;
{
  double sum=0.0;
  int N;
  int dN=20;
  int precision_OK=0;
  int j;
  fourier_int_data params;
  gsl_sum_levin_u_workspace *workspace=gsl_sum_levin_u_alloc (1+Nzeros);
  gsl_integration_workspace *gsl_wksp=gsl_integration_workspace_alloc(4000);
  gsl_function F;
  double *tmp=(double *)malloc(Nzeros*sizeof(double));
  double error;
  double *zeros;

  params.x=x;
  params.params=param;
  params.function=func;

  if (n==0){
    F.function=fourier_J0_int;
    zeros=J0zero;
  } else if (n==1){
    F.function=fourier_J1_int;
    zeros=J1zero;
  } else {
    F.function=fourier_J0_int;
    zeros=J0zero;
  }
 
  F.params=&params;
  j=Ni;
  N=dN;
  do {
    for(;j<N+Ni;j++){
      tmp[j-Ni]=integrate(zeros[j]/x,zeros[j+1]/x,F,gsl_wksp);
    }
    gsl_sum_levin_u_accel(tmp,N,workspace,&sum,&error);
    if (workspace->terms_used<N) {
      precision_OK=1;
    } else {
      if (fabs(error/sum)<=epsilon) {
	precision_OK=1;
      } else {
	N+=dN;
      }
    }
  } while((precision_OK==0)&&(N<=Nzeros));
  
  if (workspace->terms_used>=Nzeros){
    if (fabs(error/sum)>epsilon){
      fprintf(stderr,"Warning: terms used=%d rel_err=%.10e abs_err=%.10e\n",\
	      workspace->terms_used,fabs(error/sum),fabs(error));
    }
  }
  gsl_integration_workspace_free(gsl_wksp);
  free(tmp);
  gsl_sum_levin_u_free(workspace);
  
  return sum;
}


double fourier_j0(x,func,param)
     double x;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
{
  return fourier_j_i(x,0,func,param,0);
}


double fourier_j1(x,func,param)
     double x;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
{
  return fourier_j_i(x,0,func,param,1);
}


double fourier_j0_i(x,Ni,func,param)
     double x;
     int Ni;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
{
  return fourier_j_i(x,Ni,func,param,0);
}


double fourier_j1_i(x,Ni,func,param)
     double x;
     int Ni;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
{
  return fourier_j_i(x,Ni,func,param,1);
}



// 2-d Fourier integral of the function "func", starting at the Ni-th
// zero and ending at the Nf-th zero of the Bessel function J0

double fourier_j_if(x,Ni,Nf,func,param,n)
     double x;
     int Ni,Nf;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
     int n;
{
  double sum=0.0;
  int j;
  fourier_int_data params;
  gsl_integration_workspace *gsl_wksp=gsl_integration_workspace_alloc(4000);
  gsl_function F;
  double *zeros;
  
  params.x=x;
  params.params=param;
  params.function=func;
  
  if (n==0){
    F.function=fourier_J0_int;
    zeros=J0zero;
  } else if (n==1){
    F.function=fourier_J1_int;
    zeros=J1zero;
  } else {
    F.function=fourier_J0_int;
    zeros=J0zero;
  }

  F.params=&params;
  
  // We have to sum a finite number of terms ==> we must not use sum
  // acceleration algorithms !!
  
  for(j=Ni;j<Nf;j++){
    sum+=integrate(zeros[j]/x,zeros[j+1]/x,F,gsl_wksp);
  }
  
  gsl_integration_workspace_free(gsl_wksp);
  return sum;
}


double fourier_j0_f(x,Nf,func,param)
     double x;
     int Nf;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
{
  return fourier_j_if(x,0,Nf,func,param,0);
}


double fourier_j1_f(x,Nf,func,param)
     double x;
     int Nf;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
{
  return fourier_j_if(x,0,Nf,func,param,1);
}


double fourier_j0_if(x,Ni,Nf,func,param)
     double x;
     int Ni,Nf;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
{
  return fourier_j_if(x,Ni,Nf,func,param,0);
}


double fourier_j1_if(x,Ni,Nf,func,param)
     double x;
     int Ni,Nf;
#ifndef FORTRAN
     double (*func)(double,void*);
#else
     double (*func)(double*,void*);
#endif
     void *param;
{
  return fourier_j_if(x,Ni,Nf,func,param,1);
}


// Initialisation

void init_workspace_fourier(int N){
  Nzeros=N;
  //workspace=gsl_sum_levin_u_alloc (1+N);
}


void set_fourier_precision(double e,double e1){
  epsilon=e;
  epsilon1=e1;
}


/* This routine sets the processor register that controls the FPU mode
   -- It is set to use EXTENDED double precision, to round to the
   NEAREST floating point number, and to NOT generate an interrupt on
   floating point exceptions -- This is specific to x86 processors */

void set_fpu_state(void){
#ifdef X86
  fpu_control_t cw;
  cw=_FPU_EXTENDED+_FPU_RC_NEAREST \
    +_FPU_MASK_IM+_FPU_MASK_DM+_FPU_MASK_ZM \
    +_FPU_MASK_OM+_FPU_MASK_UM+_FPU_MASK_PM;
  _FPU_SETCW(cw);
  _FPU_GETCW(cw);
#endif
}
