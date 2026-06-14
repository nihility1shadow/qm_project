#ifndef _YYY_QHO_
#define _YYY_QHO_

#include <yyy_inlines.h>

class QHO {
private:
  int Neig;
  int npt;
  double freq;
  double mass;
  double x0;
  double   **Xmtx;
  dcomplex **Pmtx;
public:
  QHO(const int Ne, const double m, const double w) {
    freq = w;
    mass = m;
    Neig = Ne;
    Xmtx = array2d<double  >(Ne, Ne);
    Pmtx = array2d<dcomplex>(Ne, Ne);
    for(int i=0; i<Ne; i++) {
      if(i != Ne-1) {
        Xmtx[i][i+1] = sqrt(0.5*(1.0+i)/(mass*freq));
        Pmtx[i][i+1] = -I*sqrt(0.5*mass*freq*(1.0+i));
      }
      if(i) {
        Xmtx[i][i-1] = sqrt(0.5*i/(mass*freq));
        Pmtx[i][i-1] = I*sqrt(0.5*mass*freq*i);
      }
    }
  }
  ~QHO(){
    free2d(Xmtx);
    free2d(Pmtx);
  }
  void set_delx(const double delx) {x0 = delx; return;}
  void Enav(dcomplex *psi, const double Fold, const double Fscl, double *rlt) const ;
  double get_En(dcomplex *psi, double **const mixing=NULL) const ;
};

double Dalp_En(const int m, const int n, const double mss, const double w, const double delx) ;
void Enoverlap(const int Neig, dcomplex *const psi0, dcomplex *const psi1, 
    const double Fold, const double Fscl, double sigma[3]) ;

void clean_CS2En();
void prepare_CS2En(const int Neig) ;
int wfCS2En(const int Neig, const dcomplex wgt, const dcomplex alp, dcomplex *aux, dcomplex *wf) ;
void cswf(const int npt, const double xmin, const double dx, 
    const double m, const double w, const dcomplex wgt, const dcomplex alp, dcomplex *wf) ;
void wfEn2grid(const int Neig, dcomplex *const Cn, 
    const double m, const double w, const double xc,
    const int npt, const double xmin, const double dx, dcomplex *wf) ;
double ahmpes(const double x, void *v) ;

void test_CS2En() ;
double chk_mixing(const int k, const int l, const double m, const double w, const int npt, 
    const double xmin, const double dx, double **const wl, double **const wr) ;
void get_Enwf(const int Neig, const double m, const double w, const double xc,
    const int npt, const double xmin, const double dx, double **wf) ;

#endif
