#include <gsl/gsl_sf_lambert.h>

#include <sys/times.h>
#include <sys/time.h>

#ifdef _YYY_MPI_
#include <mpi.h>
//#include "./yyy_mpi.h"
#endif

#include <matrix.h>

#include "./ahm.h"
#include "./qj.h"
#include "./qmfft.h"
#include "./qho.h"

#define _YYY_EXACT_U1_

static void H1psi(const int dim, const double t, const double Fa, const double Fdt,
                void *pm, dcomplex *const yin, dcomplex *dydt) {
  AHM *model = (AHM *)pm;
  model->exwf(dim, Fa, Fdt, yin, dydt);
  return;
}

/* 
 * calculate |phi> = -I H1 |psi>
 * H1 = \sum_k cpl[k] (\hat{c}_k^\dagger \hat{d} + \hat{d}^\dagger\hat{c}_k)
 */
void AHM::exwf(const int dim, const double Fa, const double Fdt, 
    dcomplex *const psi, dcomplex *phi) const {
  if(Nel>1) {
    MBexwf(dim, Fa, Fdt, psi, phi);
    return;
  }

  dcomplex *qv = array1d<dcomplex>(Norb);
  for(int j=0; j<Norb; j++) qv[j] = -I*Fdt*cpl[j];

  for(int k=0; k<npt; k++) {
    phi[k] *= Fa;
    for(int j=1; j<Norb; j++) {
      int idx =j*npt+k; 
      phi[k]   += qv[j]*psi[idx];
      phi[idx]  = phi[idx]*Fa + qv[j]*psi[k];
    }
  }

  free1d(qv);
  return;
}

/* 
 * Propgation Exp(-I*H1*dt) with
 * H1 = \sum_k cpl[k] (\hat{c}_k^\dagger \hat{d} + \hat{d}^\dagger\hat{c}_k)
 */
void AHM::U1(const double dt, dcomplex **psi) const {
  dcomplex *kap = array1d<dcomplex>(Nhs*npt);
  Clsrk8(*psi, kap, Nhs*npt, 0, dt, (void *)this, H1psi);
  free1d(kap);
  return;
}

/*
 * Quamtum mechanical propagation of the Anderson-Holstein model
 * wf[n][k]: n = 0, ..., Norb-1
 *           k = 0, ..., npt-1;
 *          wf[0] is for the excited state |e>
 *          wf[k] (k>0) is for the excited state |g>|1>_{k-1}
 *          Wavefunctions are with grid representation
 */
void AHM::qm(const int nstep, const double dt, const dcomplex alp0) const {
  const int nstate = Nhs; 
  const double dx  = (xmax-xmin)/npt;
  dcomplex **wf    = array2d<dcomplex>(nstate, npt);
  double   **pes   = array2d<double>(nstate, npt);

  struct tms st_cpu, en_cpu, pr_cpu;
  clock_t st_time = times(&st_cpu);

#ifdef _YYY_EXACT_U1_
  const int use_exact_u1 = (Nhs == Norb);
  dcomplex *psi = use_exact_u1 ? array1d<dcomplex>(nstate) : 0,
           *phi = use_exact_u1 ? array1d<dcomplex>(nstate) : 0;
  CMatrix U1t = calc_U1t(dt);
#endif

  for(int j=0; j<npt; j++) {
    double x  = xmin + j*dx,
           hp = 0.5*mass*freq*freq*x*x;
    for(int n=0;   n<Nex;    n++) pes[n][j] = 0.5*mass*freq*freq*(x-delx)*(x-delx) + Eocc[n];
    for(int n=Nex; n<nstate; n++) pes[n][j] = hp + Eocc[n];
  }

  //initial state
  cswf(npt, xmin, dx, mass, freq, 1.0, alp0, wf[0]);
 
  char fnm[256];
  sprintf(fnm, "ahm-qm-s%d-n%d.dat", Norb, Nel);
  FILE *FL = fopen(fnm, "w");

  fprintf(FL, "#discretizing the bath\n");
  for(int n=0; n<Norb; n++) fprintf(FL, "#%6d %1.16e %1.16e\n", n, cpl[n], En[n]);
  wfanalysis(FL, 0.0, wf);

  QMFFT1D sys(npt, xmin, xmax);
  sys.set_mass(mass);


  clock_t pr_time = times(&pr_cpu);
  for(int t=0; t<nstep; t++) {
    for(int n=0; n<nstate; n++) sys.propagator(0.5*dt, pes[n], wf[n]);
#ifdef _YYY_EXACT_U1_
    if(use_exact_u1) {
      for(int j=0; j<npt; j++) {
        for(int n=0; n<nstate; n++) psi[n] = wf[n][j];
        U1t.vectorply(psi, phi);
        for(int n=0; n<nstate; n++) wf[n][j] = phi[n];
      }
    } else {
      U1(dt, wf);
    }
#else
    U1(dt, wf);
#endif
    for(int n=0; n<nstate; n++) sys.propagator(0.5*dt, pes[n], wf[n]);
    wfanalysis(FL, (t+1)*dt, wf);
  }
  clock_t en_time = times(&en_cpu);
  print_timing(st_time, pr_time, en_time,
               st_cpu,  pr_cpu,  en_cpu, FL);
  fclose(FL);


  sprintf(fnm, "ahmwf-qm-s%d.dat", Norb);
  FL = fopen(fnm, "w");
  //dcomplex alp0 = sqrt(mass*freq/2)*(xinit+(pinit/(mass*freq))*I);
  //double p0 =  real(alp0)*mass*freq/sqrt(mass*freq/2), 
  //       q0 = real(alp0)/sqrt(mass*freq/2);
  //dcomplex ph = exp(-0.5*p0*q0*I);
  for(int j=0; j<npt; j++) {
    double x = xmin + j*dx;
    fprintf(FL, "%1.8f", x);
    for(int n=0; n<nstate; n++) fprintf(FL, " %+1.16e %+1.16e", real(wf[n][j]), imag(wf[n][j]));
    fprintf(FL, "\n");
  }
  fclose(FL);

  free2d(wf );
  free2d(pes);

#ifdef _YYY_EXACT_U1_
  if(use_exact_u1) {
    free1d(psi);
    free1d(phi);
  }
#endif
  return;
}

/*
 * MPI version of AHM::qm
 * Quamtum mechanical propagation of the Anderson-Holstein model
 * wf[n][k]: n = 0, ..., Norb-1
 *           k = 0, ..., npt-1;
 *          wf[0] is for the excited state |e>
 *          wf[k] (k>0) is for the excited state |g>|1>_{k-1}
 *          Wavefunctions are with grid representation
 */
void AHM::qm2(const int nstep, const double dt, const dcomplex alp0) const {
  const int nstate = Nhs; 
  const double dx  = (xmax-xmin)/npt;
  dcomplex **wf    = array2d<dcomplex>(nstate, npt);
  double   **pes   = array2d<double>(nstate, npt);

  for(int j=0; j<npt; j++) {
    double x  = xmin + j*dx,
           hp = 0.5*mass*freq*freq*x*x;
    for(int n=0;   n<Nex;    n++) pes[n][j] = 0.5*mass*freq*freq*(x-delx)*(x-delx) + Eocc[n];
    for(int n=Nex; n<nstate; n++) pes[n][j] = hp + Eocc[n];
  }

  //initial state
  cswf(npt, xmin, dx, mass, freq, 1.0, alp0, wf[0]);
 
  char fnm[256];
  sprintf(fnm, "ahm-qm-s%d-n%d.dat", Norb, Nel);
  FILE *FL = fopen(fnm, "w");

  fprintf(FL, "#discretizing the bath\n");
  for(int n=0; n<Norb; n++) fprintf(FL, "#%6d %1.16e %1.16e\n", n, cpl[n], En[n]);
  wfanalysis(FL, 0.0, wf);

  QMFFT1D sys(npt, xmin, xmax);
  sys.set_mass(mass);

  for(int t=0; t<nstep; t++) {
    for(int n=0; n<nstate; n++) sys.propagator(0.5*dt, pes[n], wf[n]);
    U1(dt, wf);
    for(int n=0; n<nstate; n++) sys.propagator(0.5*dt, pes[n], wf[n]);
    wfanalysis(FL, (t+1)*dt, wf);
  }
  fclose(FL);

  sprintf(fnm, "ahmwf-qm-s%d.dat", Norb);
  FL = fopen(fnm, "w");
  //dcomplex alp0 = sqrt(mass*freq/2)*(xinit+(pinit/(mass*freq))*I);
  //double p0 =  real(alp0)*mass*freq/sqrt(mass*freq/2), 
  //       q0 = real(alp0)/sqrt(mass*freq/2);
  //dcomplex ph = exp(-0.5*p0*q0*I);
  for(int j=0; j<npt; j++) {
    double x = xmin + j*dx;
    fprintf(FL, "%1.8f", x);
    for(int n=0; n<nstate; n++) fprintf(FL, " %+1.16e %+1.16e", real(wf[n][j]), imag(wf[n][j]));
    fprintf(FL, "\n");
  }
  fclose(FL);

  free2d(wf );
  free2d(pes);

  return;
}

/*
 * rlt[0] : average number of electrons
 * rlt[1] : average position
 * rlt[2] : average vibration energy
 * rlt[3+n] (0<n<Norb) : average occupation in the n-th orbital
 *
 * the argument excited indicates the state is excited or not.
 * Though this argument is redundant because it can be directly deduced from state, 
 * it spares the caculation since it is already known when this function is called.
 */
void AHM::csproj(const dcomplex wgt_for, const dcomplex alp_for, 
    const dcomplex wgt_bck, const dcomplex alp_bck,
    const int excited, const set<int> &state, double *rlt) const {
  static const double sqrt1_hmw = sqrt(2.0/(mass*freq)), 
                      sqrt2mw   = sqrt(2.0*mass*freq),
                      inv_2m    = 1/(2*mass),
                      hw        = 0.5*freq,
                      m2        = mass*mass,
                      w2        = freq*freq,
                      dmw       = 2*mass*freq,
                      hmw2      = 0.5*mass*freq*freq,
                      inv_4mw   = 1/(4*mass*freq),
                      xshift[2] = {0, delx};

  double pf = imag(alp_for)*sqrt2mw  ,
         xf = real(alp_for)*sqrt1_hmw,
         pb = imag(alp_bck)*sqrt2mw  ,
         xb = real(alp_bck)*sqrt1_hmw;
  dcomplex ovlp = exp(-inv_4mw*((pf-pb)*(pf-pb) + m2*w2*(xf-xb)*(xf-xb) - (dmw*(pf*xb-pb*xf))*I)),
           prb  = wgt_for*conj(wgt_bck)*ovlp;
  double prb_rl = real(prb);
  for(int j : state) {
    rlt[3+j] += prb_rl;
    rlt[0]   += prb_rl;
  }

  dcomplex xav = 0.5*(xf+xb) + (2*inv_4mw*(pf-pb))*I,
           pav = 0.5*(pf+pb) + (0.25*dmw*(xb-xf))*I;
  rlt[1] += real(prb*xav);
  rlt[2] += real(prb*(hw + pav*pav*inv_2m + hmw2*(xav-xshift[excited])*(xav-xshift[excited])));
  return;
}

void AHM::wfanalysis(FILE *FL, const double t, dcomplex **const wf) const {
  double rlt[3], pm[3], sum, xav, ven, 
         *prb=array1d<double>(Norb);
  pm[0] = freq;
  pm[1] = mass;
  pm[2] = delx;

  QMFFT1D sys(npt, xmin, xmax);
  sys.set_mass(mass);

#if 0
  if(Nel==1) {
    qav(xmin, xmax, npt, wf[0], 0.0, 1.0, &rlt[0]);
    rlt[2] = sys.get_En(mass, (void *)&pm[0], wf[0], &ahmpes);
    fprintf(FL, "%12.8f %+1.16e %+1.16e %+1.16e\n", t, 
      rlt[0], rlt[1]/rlt[0], rlt[2]/rlt[0]);
  } else 
#endif
  {
    xav = 0.0;
    ven = 0.0;
    sum = 0.0;
    bzero(prb, sizeof(double)*Norb);
    for(int n=0; n<Nhs; n++) {
      qav(xmin, xmax, npt, wf[n], 0.0, 1.0, &rlt[0]);
      if(n<Nex) pm[2] = delx;
      else      pm[2] = 0.0;
      rlt[2] = sys.get_En(mass, (void *)&pm[0], wf[n], &ahmpes);
      sum += rlt[0];
      xav += rlt[1];
      ven += rlt[2];
      for(int k=0; k<Nel; k++) prb[occ[n][k]] += rlt[0];
    }
    fprintf(FL, "%12.8f %+1.16e %+1.16e %+1.16e", t, sum, xav, ven);
    for(int k=0; k<Norb; k++) fprintf(FL, " %+1.16e", prb[k]);
    fprintf(FL, "\n");
  }

  return;
}

void AHM::discretize(const int Nstate, const double eta, const double wc){
  Norb = Nstate; 
  cpl = array1d<double>(Norb);
  En  = array1d<double>(Norb);
  double c = wc*sqrt(eta/(Norb-1)), tm0;

  En[1]  = wc/(1.32*sqrt(Norb-1.0));
  cpl[1] = c;
  tm0 = 1+En[1]/wc;
  tm0 = tm0*exp(-tm0);
  for(int n=1; n<Norb-1; n++) {
    double tmp = tm0 - n*(1-En[1]*En[1]/(wc*wc))/((Norb-2)*M_E);
    cpl[1+n] = c;
    En[1+n]  = (-gsl_sf_lambert_W0(-tmp)-1)*wc;
    if(En[1+n]<0) {
      En[1+n]  = (-gsl_sf_lambert_Wm1(-tmp)-1)*wc;
    }
    if(En[1+n]<0) {
      printf("failed to discretize the bath.\n");
      abort();
    }
  }
  
  char fnm[256];
  sprintf(fnm, "bath-dis-%d.dat", Norb);
  FILE *FL = fopen(fnm, "w");
  double sum = 0.0;
  fprintf(FL, "#discretizing the bath eta = %g and wc = %g\n", eta, wc);
  for(int n=1; n<Norb; n++) {
    double Iw = eta*wc*wc*(1-(1+En[n]/wc)*exp(-En[n]/wc));
    fprintf(FL, "%6d %1.16e %1.16e %1.16e\n", n, c, En[n], Iw);
    sum += wc/En[n];
  }
  fprintf(FL, "#The renormalization energy is %g eta*wc, the target is eta*wc\n", sum/(Norb-1));
  fclose(FL);

  return;
}

void AHM::diseven2(const int Nstate, const double eta, const double wc){
  Norb = Nstate; 
  cpl = array1d<double>(Nstate);
  En  = array1d<double>(Nstate);
  double c = eta/sqrt(Norb-1);

  for(int n=0; n<Norb-1; n++) {
    cpl[n+1] = c;
    En[n+1]  = -(Norb-2.0-n)*wc/(Norb-2);
  }

  return;
}

void AHM::diseven(const int Nstate, const double eta, const double wc){
  Norb = Nstate; 
  cpl = array1d<double>(Nstate);
  En  = array1d<double>(Nstate);

  if(Norb < 3) {
    printf("failed to build bath: Norb must be at least 3.\n");
    abort();
  }

  const double eV_to_au = 1.0/27.211386245988;
  const double EF = -4.5 * eV_to_au;
  const double Eb = -10.0 * eV_to_au;

  const double E0 = EF - 0.5*wc;
  const double E1 = EF + 0.5*wc;

  if(E0 <= Eb || E1 <= E0) {
    printf("failed to build bath.\n");
    printf("Need Eb < E0 < E1, but Eb = %g, E0 = %g, E1 = %g.\n", Eb, E0, E1);
    abort();
  }

  const double A = E0 - Eb;
  const double B = E1 - Eb;
  const double M = Norb - 1.0;
  const double r = pow(B/A, 2.0/3.0);

  if(r <= 1.0) {
    printf("failed to build bath: invalid endpoint ratio r = %g.\n", r);
    abort();
  }

  const double joff = (M - r)/(r - 1.0);

  if(joff <= -1.0) {
    printf("failed to build bath: invalid joff = %g.\n", joff);
    abort();
  }

  const double k = A/pow(1.0 + joff, 1.5);
  const double c = sqrt(eta/(Norb - 1.0));

  En[0]  = 0.0;
  cpl[0] = 0.0;

  for(int n=0; n<Norb-1; n++) {
    const double j = n + 1.0;
    En[n+1]  = Eb + pow(j + joff, 1.5)*k;
    cpl[n+1] = c;
  }

  return;
}

void AHM::disrandom(const int Nstate, const double eta, const double wc) {
  Norb = Nstate; 
  cpl = array1d<double>(Nstate);
  En  = array1d<double>(Nstate);

  for(int n=0; n<Norb; n++) {
    En[1+n]  = (n+0.5)*wc/(Norb-1);
    cpl[1+n] = sqrt(eta*En[1+n]/(Norb-1))*exp(-0.5*En[1+n]/wc);
  }

  return;
}
