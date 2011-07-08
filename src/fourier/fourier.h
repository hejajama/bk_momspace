// Last modified on June 06th, 2006 -- F. Gelis


extern int Nzeros;
extern double *workspace;


#ifndef FORTRAN
double fourier_j0(double x,double (*func)(double,void*),void* param);
double fourier_j0_i(double x,int Ni,double (*func)(double,void*),void *param);
double fourier_j0_f(double x,int Nf,double (*func)(double,void*),void *param);
double fourier_j0_if(double x,int Ni,int Nf,double (*func)(double,void*),void *param);
double fourier_j1(double x,double (*func)(double,void*),void* param);
double fourier_j1_i(double x,int Ni,double (*func)(double,void*),void *param);
double fourier_j1_f(double x,int Nf,double (*func)(double,void*),void *param);
double fourier_j1_if(double x,int Ni,int Nf,double (*func)(double,void*),void *param);
#else
double fourier_j0(double x,double (*func)(double*,void*),void* param);
double fourier_j0_i(double x,int Ni,double (*func)(double*,void*),void *param);
double fourier_j0_f(double x,int Nf,double (*func)(double*,void*),void *param);
double fourier_j0_if(double x,int Ni,int Nf,double (*func)(double*,void*),void *param);
double fourier_j1(double x,double (*func)(double*,void*),void* param);
double fourier_j1_i(double x,int Ni,double (*func)(double*,void*),void *param);
double fourier_j1_f(double x,int Nf,double (*func)(double*,void*),void *param);
double fourier_j1_if(double x,int Ni,int Nf,double (*func)(double*,void*),void *param);
#endif

void init_workspace_fourier(int N);
void set_fpu_state(void);
void set_fourier_precision(double e,double e1);

