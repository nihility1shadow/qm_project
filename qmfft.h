#ifndef _YYY_QM_FFTW_
#define _YYY_QM_FFTW_

#include <fftw3.h>
#include <yyy_inlines.h>

typedef double RF2p(const double, void *);
typedef dcomplex CF2p(const double, void *);

class QMFFT1D {
  private: 
    int npt;
    dcomplex *psi; 
    dcomplex *kap; 
    fftw_plan fftpln_for; 
    fftw_plan fftpln_bck;
    double xmin;
    double xmax;
    double mass;
  public:
    QMFFT1D(){
      npt = 0;
      psi = NULL;
      kap = NULL;
    }
    ~QMFFT1D(){
      free1d(psi);
      free1d(kap);
      fftw_destroy_plan(fftpln_for);
      fftw_destroy_plan(fftpln_bck);
    }
    QMFFT1D(const int xstep, const double x0, const double x1) {
      npt = xstep;
      xmin = x0;
      xmax = x1;
      psi = array1d<dcomplex>(npt);
      kap = array1d<dcomplex>(npt);
      fftpln_for = fftw_plan_dft_1d(npt, 
                  (fftw_complex *)psi, (fftw_complex *)kap, 
                  FFTW_FORWARD,  FFTW_MEASURE);
      fftpln_bck = fftw_plan_dft_1d(npt, 
                  (fftw_complex *)kap, (fftw_complex *)psi, 
                  FFTW_BACKWARD, FFTW_MEASURE);     
    }
    void set_mass(const double m) {mass = m;}
    void set_psi(dcomplex *const data){
      memcpy(psi, data, sizeof(dcomplex)*npt);
      return;
    }
    void get_psi(dcomplex *data) const{
      memcpy(data, psi, sizeof(dcomplex)*npt);
      return;
    }
    void propagator(const double dt, const double *v, dcomplex *wf);
    void pnpsi(const int n, dcomplex *pnpsi);
    double get_En(const double mass, void *pm, dcomplex *psi, double (*v)(const double, void *));
};

class Sys2e {
  public:
    int npt;
    double *Eg;
    double *Ee;
    double *Vna;   //nonadiabatic term <g| H_{nonadiabatic} | e>
    dcomplex *Veg;
    dcomplex *Vge;
    double xmin;
    double xmax;
    double mass;
    QMFFT1D qft;
    Sys2e(const int n, const double x0, const double x1, const double m0) : qft(n, x0, x1){
      npt  = n;
      xmin = x0;
      xmax = x1;
      mass = m0;
      Eg   = array1d<double>(npt);
      Ee   = array1d<double>(npt);
      Vna  = array1d<double>(npt);
      Vge  = array1d<dcomplex>(npt);
      Veg  = array1d<dcomplex>(npt);
    }
    ~Sys2e() {
      free1d(Eg);
      free1d(Ee);
      free1d(Vna);
      free1d(Veg);
      free1d(Vge);
    }
    void set_Eg(double *const En){
      memcpy(Eg, En, sizeof(double)*npt);
      return;
    }
    void set_Eg(void *pm, RF2p *f){
      double dx = (xmax-xmin)/npt;
      for(int j=0; j<npt; j++) {
        double x = xmin + j*dx;
        Eg[j] = (*f)(x, pm);
      }
      return;
    }
    void set_Ee(void *pm, RF2p *f){
      double dx = (xmax-xmin)/npt;
      for(int j=0; j<npt; j++) {
        double x = xmin + j*dx;
        Ee[j] = (*f)(x, pm);
      }
      return;
    }
    void set_Ee(double *const En){
      memcpy(Ee, En, sizeof(double)*npt);
      return;
    }
    void set_Vna(void *pm, RF2p *f){
      double dx = (xmax-xmin)/npt;
      for(int j=0; j<npt; j++) {
        double x = xmin + j*dx;
        Vna[j] = (*f)(x, pm);
      }
      return;
    }
    void set_Vna(double *const En){
      memcpy(Vna, En, sizeof(double)*npt);
      return;
    }
    void set_Veg(void *pm, CF2p *f){
      double dx = (xmax-xmin)/npt;
      for(int j=0; j<npt; j++) {
        double x = xmin + j*dx;
        Veg[j] = (*f)(x, pm);
      }
      return;
    }
    void set_Veg(dcomplex *const En){
      memcpy(Veg, En, sizeof(dcomplex)*npt);
      return;
    }
    void set_Vge(void *pm, CF2p *f){
      double dx = (xmax-xmin)/npt;
      for(int j=0; j<npt; j++) {
        double x = xmin + j*dx;
        Vge[j] = (*f)(x, pm);
      }
      return;
    }
    void set_Vge(dcomplex *const En){
      memcpy(Vge, En, sizeof(dcomplex)*npt);
      return;
    }
    void hermicity();
};

#endif
