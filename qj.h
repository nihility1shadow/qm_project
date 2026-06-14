#ifndef _YYY_POISSON_REP_H_
#define _YYY_POISSON_REP_H_

#include <fftw3.h>
#include <yyy_inlines.h>
#include <myid.h>

#define _MODIFYING_TULLY1_ 1

typedef double _Pes(const double);
typedef double _RF2p(const double, const double);
typedef double _RF4p(const double x, const double h, const double x0, const double width);

int test_adiabatic(int argc, char **argv) ;
int adiabatic_hk(int argc, char **argv) ;
int  diabatic(int argc, char **argv) ;
int  repeigen(int argc, char **argv) ;
void check_df(const int npt, const double xmin, const double dx, 
    double (*f)(const double), double (*df)(const double), char *);
void check_df(const int npt, const double xmin, const double dx, const double pm,
    double (*f)(const double, const double), double (*df)(const double, const double), char *);
void check_df(const int npt, const double xmin, const double x, void *pm,
    double (*f)(const double, void *), double (*df)(const double, void *), char *);

dcomplex init_gaussian(const double x, const double x0, const double k0, const double sigma) ;
double gaussian_shape(const double x, const double height, const double x0, const double width) ;

void qav(const double xmin, const double xmax, const int xstep, 
			dcomplex *psi, const double Fold, const double Fscl, double *rlt) ;
void correlation(const double xmin, const double xmax, const int xstep, 
    dcomplex *psi0, dcomplex *psi1, 
    const double Fold, const double Fscl, double sigma[3]) ;

void PJ_tls(const int ntraj, const double omega, const double kappa, 
    const double rtfct, const int nstep, const double dt);
void BP_tls(const int ntraj, const double omega, const double kappa, 
    const double rtfct, const int nstep, const double dt);
void PJ_tls_rr(const int ntraj, const double omega, const double kappa, 
    const int nstep, const double dt) ;

void BP_const(const int ntraj, const double rtfct, const int nstep, const double dt) ;

void QJpoisson(const int ntraj, double *const lambda, 
      const int nstep, const double dt, 
      const int nwf,   dcomplex **wf1, dcomplex **wf2, double **sigma,
      const double xmin, const double xmax, const int xstep, 
			const double mass, double *const v1, double *const v2, double *const vc,
      dcomplex *psi1, dcomplex *kap1, 
      fftw_plan fftpln_for1, fftw_plan fftpln_bck1,
      dcomplex *psi2, dcomplex *kap2, 
      fftw_plan fftpln_for2, fftw_plan fftpln_bck2) ;

void NHpoisson(const int ntraj, double *const lambda, 
      const int nstep, const double dt, 
      const int nwf,   dcomplex **wf1, dcomplex **wf2, double **sigma,
      const double xmin, const double xmax, const int xstep, 
			const double mass, double *const v1, double *const v2, double *const vc, 
      dcomplex *psi1, dcomplex *kap1, 
      fftw_plan fftpln_for1, fftw_plan fftpln_bck1,
      dcomplex *psi2, dcomplex *kap2, 
      fftw_plan fftpln_for2, fftw_plan fftpln_bck2);

void NHpoisson_norepeat(const int ntraj, double *const lambda, 
      const int nstep, const double dt, 
      const int nwf,   dcomplex **wf1, dcomplex **wf2, double **sigma,
      const double xmin, const double xmax, const int xstep, 
			const double mass, double *const v1, double *const v2, double *const vc, 
      dcomplex *psi1, dcomplex *kap1, 
      fftw_plan fftpln_for1, fftw_plan fftpln_bck1,
      dcomplex *psi2, dcomplex *kap2, 
      fftw_plan fftpln_for2, fftw_plan fftpln_bck2);

void QJqm(const double xmin, const double xmax, const int xstep, 
			const double t, const double dt, const double mass, 
      double *const v1, double *const v2, double *const vc,
      dcomplex *psi1, dcomplex *kap1, 
      fftw_plan fftpln_for1, fftw_plan fftpln_bck1,
      dcomplex *psi2, dcomplex *kap2, 
      fftw_plan fftpln_for2, fftw_plan fftpln_bck2);

void RabiPoisson(const int ntraj, const int nstep, const double dt, 
    const int npt, const double xmin, const double xmax,
    const double delx, const double delE, const double lambda,
    const double xinit);
void RabiNHPP(const int ntraj, const int nstep, const double dt, 
    const int npt, const double xmin, const double xmax,
    const double delx, const double delE, const double lambda,
    const double xinit) ;

void adiapoisson(const int ntraj, const int nstep, const double dt, 
      dcomplex **wf1, dcomplex **wf2, double **sigma,
      const double xmin, const double xmax, const int xstep, 
			const double mass, double *const v1, double *const v2, 
      double *const kij, double *const dij, 
      dcomplex *psi1, dcomplex *psi2); 

double Fzero1p(const double x) ;
double Fzero2p(const double x, const double kappa) ;
double Hint_const(const double x, const double kappa);
double coh_constint(const double q0, const double p0, const double gmm, const double kappa);
double coh_constint(const double q0, const double p0, const double gmm, const double kappa) ;
double coh_tully1int(const double q, const double p, const double gmm, const double kappa) ;
double coh_tully2int(const double q, const double p, const double gmm, const double kappa) ;
double coh_tully3int(const double q, const double p, const double gmm, const double kappa) ;
double coh_tanhint(const double q, const double p, const double gmm, const double kappa) ;


double     tully1up(const double x) ;
double der_tully1up(const double x) ;
double  d2_tully1up(const double x) ;
double  d3_tully1up(const double x) ;
double     tully1down(const double x) ;
double der_tully1down(const double x) ;
double  d2_tully1down(const double x) ;
double  d3_tully1down(const double x) ;
double     tully1int(const double x, const double pm) ;
double der_tully1int(const double x, const double pm) ;
double  d2_tully1int(const double x, const double pm) ;
double  d3_tully1int(const double x, const double pm) ;
void dVall_tully1up(const double x, double *f); 
void dVall_tully1down(const double x, double *f); 
void dHall_tully1(const double x, const double pm, double *f); 

double     tully2up(const double x) ;
double der_tully2up(const double x) ;
double  d2_tully2up(const double x) ;
double  d3_tully2up(const double x) ;
double     tully2down(const double x) ;
double der_tully2down(const double x) ;
double  d2_tully2down(const double x) ;
double  d3_tully2down(const double x) ;
double     tully2int(const double x, const double pm) ;
double der_tully2int(const double x, const double pm) ;
double  d2_tully2int(const double x, const double pm) ;
double  d3_tully2int(const double x, const double pm) ;
void dVall_tully2up(const double x, double *f); 
void dVall_tully2down(const double x, double *f); 
void dHall_tully2(const double x, const double pm, double *f); 

double     tully3up(const double x) ;
double der_tully3up(const double x) ;
double  d2_tully3up(const double x) ;
double  d3_tully3up(const double x) ;
double     tully3down(const double x) ;
double der_tully3down(const double x) ;
double  d2_tully3down(const double x) ;
double  d3_tully3down(const double x) ;
double     tully3int(const double x, const double pm) ;
double der_tully3int(const double x, const double pm) ;
double  d2_tully3int(const double x, const double pm) ;
double  d3_tully3int(const double x, const double pm) ;
void dVall_tully3up(const double x, double *f); 
void dVall_tully3down(const double x, double *f); 
void dHall_tully3(const double x, const double pm, double *f); 

double     LZdown(const double x) ;
double der_LZdown(const double x) ;
double  d2_LZdown(const double x) ;
double  d3_LZdown(const double x) ;
double     LZup(const double x) ;
double der_LZup(const double x) ;
double  d2_LZup(const double x) ;
double  d3_LZup(const double x) ;
void dVall_LZup(const double x, double *f); 
void dVall_LZdown(const double x, double *f); 
void dHall_LZ(const double x, const double pm, double *f); 

double     tanhup(const double x) ;
double der_tanhup(const double x) ;
double  d2_tanhup(const double x) ;
double  d3_tanhup(const double x) ;
double     tanhdown(const double x) ;
double der_tanhdown(const double x) ;
double  d2_tanhdown(const double x) ;
double  d3_tanhdown(const double x) ;
double     tanhint(const double x, const double pm) ;
double der_tanhint(const double x, const double pm) ;
double  d2_tanhint(const double x, const double pm) ;
double  d3_tanhint(const double x, const double pm) ;
void dVall_tanhup(const double x, double *f); 
void dVall_tanhdown(const double x, double *f); 
void dHall_tanh(const double x, const double pm, double *f); 

double     rabiup(const double x) ;
double der_rabiup(const double x) ;
double  d2_rabiup(const double x) ;
double  d3_rabiup(const double x) ;
double     rabidown(const double x) ;
double der_rabidown(const double x) ;
double  d2_rabidown(const double x) ;
double  d3_rabidown(const double x) ;
void dVall_rabiup(const double x, double *f); 
void dVall_rabidown(const double x, double *f); 
void dHall_rabi(const double x, const double pm, double *f); 

double     AHMup(const double x) ;
double der_AHMup(const double x) ;
double  d2_AHMup(const double x) ;
double  d3_AHMup(const double x) ;
double     AHMdown(const double x) ;
double der_AHMdown(const double x) ;
double  d2_AHMdown(const double x) ;
double  d3_AHMdown(const double x) ;
void dVall_AHM(const double x, double *f); 
void dHall_AHM(const double x, const double pm, double *f); 

extern double _LU_eps;
double     LUZHup(const double x) ;
double der_LUZHup(const double x) ;
double  d2_LUZHup(const double x) ;
double  d3_LUZHup(const double x) ;
double     LUZHdown(const double x) ;
double der_LUZHdown(const double x) ;
double  d2_LUZHdown(const double x) ;
double  d3_LUZHdown(const double x) ;
double     LUZHint(const double x, const double pm) ;
double der_LUZHint(const double x, const double pm) ;
double  d2_LUZHint(const double x, const double pm) ;
double  d3_LUZHint(const double x, const double pm) ;
double coh_LUZHint(const double q, const double p, const double gmm, const double kappa) ;
void dVall_LUZHup(const double x, double *f); 
void dVall_LUZHdown(const double x, double *f); 
void dHall_LUZH(const double x, const double pm, double *f); 

void DpsiDt(const int dim, const double t, const double a, const double b,
            void *pm, dcomplex *const yin, dcomplex *dydt) ;

void print_timing(
    const clock_t &st_time, 
    const clock_t &pr_time, 
    const clock_t &en_time,
    const struct tms &st_cpu,
    const struct tms &pr_cpu,
    const struct tms &en_cpu,
    FILE *const logFL);

void Clsrk8(
       dcomplex *y, 
       dcomplex *dy, 
       const int n, 
       const double x, 
       const double h, 
       void *pm, 
       void (*derivs)(const int dim,const double t, const double a, const double b, 
         void *pm, dcomplex *const yinit, dcomplex *dydt));

int fghpm(const double zmu, const int nwrite, const int Nr, const double rmin, const double rmax, 
		double *eigval, double **eigvec, double (*pes1d)(double)) ;

extern "C"
void qm1d_tdpes(const double xmin, const double xmax, const int xstep, 
			const double t, const double dt, dcomplex *psi, dcomplex *kap, 
			fftw_plan fftpln_for, fftw_plan fftpln_bck,
			const double mass, const double *v);

extern "C"
void pnpsi(const int n, const double xmin, const double xmax, const int xstep,
      dcomplex *psi, dcomplex *kap, fftw_plan fftpln_for, fftw_plan fftpln_bck) ;
#endif
