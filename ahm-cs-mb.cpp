
#include <gsl/gsl_sf_lambert.h>
#ifdef _YYY_MPI_
#include <mpi.h>
//#include "./yyy_mpi.h"
#endif

#include "./ahm.h"
#include "./qj.h"
#include "./qmfft.h"
#include "./qho.h"

#define _YYY_REALSPACE_AV_

/*
 * coherent-state approach for the Poisson method
 *
 * the coordinate representation of the coherent-state
 * <x|alpha> = (m w/(Pi hbar))^(1/4) Exp[-m w (x-xt)^2/(2 hbar) + I pt (x-xt)/hbar + I xt pt/(2 hbar)];
 * where alpha = (m w xt + I*pt)/sqrt[2 hbar m w]
 *
 * the evolution
 *
 * e^{-i H t/hbar}|alpha(0)> = Exp[-I w t/2] |alpha(t)>,
 * where (xt, pt) is trajectoy in the phase space following the classical propagation.
 *
 * In the electronic excited state, the potential is V(x) = 1/2 m w^2 (x-delx)^2,
 * xt = (x0-delx) Cos[w t] + p0/(m w) Sin[w t];
 * pt = p0 Cos[w t] - m w (x0-delx) Sin[w t];
 *
 * In other  states, the potential is V(x) = 1/2 m w^2 x^2,
 * xt = x0 Cos[w t] + p0/(m w) Sin[w t];
 * pt = p0 Cos[w t] - m w x0 Sin[w t];
 *
 */
void AHM::MBpoisson(const int ntraj, const int nstep, const double dt, 
    const dcomplex alp_init, const dcomplex amp, const int state_init) const {
  const dcomplex Iton[4] = {1, -I, -1, I};
  const int sizeofint  = sizeof(int),
            Nvac       = Norb - Nel,
            eo[2]      = {1, -1}; 
  const double dx      = -2*xmin/npt,
               fpt     = sqrt(2*mass*freq),
            sqrtNvac   = sqrt(Nvac),
            sqrt1_Nvac = 1./sqrt(Nvac),
            sqrtNorb   = sqrt(Norb),
            sqrt1_Norb = 1/sqrt(Norb),
            sqrtNel    = sqrt(Nel),
            sqrt1_Nel  = 1./sqrt(Nel),
              dalp     = mass*freq*delx/sqrt(2*mass*freq); // the shift of alpha parameter corresponding to delx
  dcomplex exphwdt     = exp(-0.5*freq*dt*I),
           expfreqdt   = exp(-freq*dt*I),
            *expEndt   = array1d<dcomplex>(Norb),
           alp0, wgt0, alp, wgt;

  for(int j=2; j<Norb; j++) if(cpl[j-1] != cpl[j]) {
    printf("At this momentum, the Poisson jump method is only valid for homogeneous cj.\n");
    abort();
  }

  for(int j=0; j<Norb; j++) expEndt[j] = exp(-I*dt*En[j]);

  int ntraj_local = ntraj, *jc = array1d<int>(nstep);
#ifdef _YYY_MPI_
  ntraj_local = ntraj/nproc;
  if(myid==0) ntraj_local += ntraj - ntraj_local*nproc;
  const double Lfct = 1.0/nproc;
#else
  const double Lfct = 1.0;
#endif

  /* the initial condition: 
   * alp0 = alp_init
   * wgt0 = 1
   */
  alp0 = alp_init;
  wgt0 = amp;
  const int nwf   = nstep/5,
            nbloc = 5;
#ifdef _YYY_REALSPACE_AV_
  dcomplex ***wf    = array3d<dcomplex>(nwf+1, Nhs, npt);
  cswf(npt, xmin, dx, mass, freq, wgt0*Lfct, alp0, wf[0][state_init]);
#else
  const int Neig = 200;
  dcomplex ***wf   = array3d<dcomplex>(nwf+1, Nhs, Neig+1),
           *aux    = array1d<dcomplex>(Neig+1);
  prepare_CS2En(Neig) ;
  wfCS2En(Neig, wgt0*Lfct, alp0, aux, wf[0][state_init]);
#endif

  const double lambda = abs(cpl[1]),
               rate   = sqrtNel*sqrtNvac*lambda*dt;
  double p0, pt, *sclf = array1d<double>(nstep+1);
  for(int j=0; j<=nstep; j++) sclf[j] = exp(rate*j)/ntraj;

  int idx, istate,
      *state = array1d<int>(Nel),      // the occupied orbitals
      *vac   = array1d<int>(Nvac),     // the vacant orbitals
      *vecn  = array1d<int>(Norb);     // an auxiliary vector
  int excited0   = occ[state_init][0] == 0 ? 1 : 0;
  for(int n=0; n<ntraj_local; n++) {
    alp    = alp0;
    wgt    = wgt0;
    istate = state_init;
    for(int j=0; j<Nel; j++) state[j] = occ[istate][j];
    complement(Norb, vecn, Nel, state, vac);

    int excited = excited0, //the index show the molecule is excited (1) or not (0)
        nj      = 0,        //number of jumps;
        sign, found;
#ifdef _TRACE_STATE_
    printf("traj: %6d\n", n);
#endif
    for(int j=1; j<=nstep; j++) {
      //printf("    time: %4d", j);
      // jump seperatedly, switch the state between 0 and idx
      if(drand48() < rate) {
        //idx in this block indicate the orbital involved in the quantum jump
#ifdef _TRACE_STATE_
        printf("    "); for(int l=0; l<Nel; l++) printf("%d", state[l]);
#endif
        if(excited) {
          // if the molecule is excited, the electron jumps to an vacant orbital
          idx = drand48()*Nvac;
          idx = vac[idx];
          state[0] = idx;
          qsort(state, Nel, sizeofint, intcmp);
          for(int l=0; l<Nel; l++) if(state[l] == idx) {
            sign = eo[l%2];
            break;
          }
          wgt *= sign*sqrtNvac*sqrt1_Nel;
#ifdef _TRACE_STATE_
          printf(" : sign = %+d, switch 0 -> %4d : ", sign, idx);
#endif
        } else {
          // if the molecule is not excited, the electron jumps from an occupied orbital to the molecule
          idx  = drand48()*Nel;
          sign = eo[idx%2];
          wgt *= sign*sqrtNel*sqrt1_Nvac;
#ifdef _TRACE_STATE_
          printf(" : sign = %+d, switch 0 <- %4d : ", sign, state[idx]);
#endif
          state[idx] = 0;
          qsort(state, Nel, sizeofint, intcmp);
        }
        nj++; 
        excited = 1 - excited;
        complement(Norb, vecn, Nel, state, vac);
        istate = binary_search3(state, *occ, Nhs, Nel*sizeof(int), &found, _mycmp3);
#ifdef _TRACE_STATE_
        for(int l=0; l<Nel; l++) printf("%d", state[l]); printf("\n"); 
#endif
      }
      /*
       * when nj is even, the electron is at state 0, 
       *     only alp[0] and wgt[0] are meaningfull, 
       *     and it is illegal to reference alp[k] and wgt[k] for k>0
       *
       * when nj is odd, the electron is at state n (n>0), 
       *     only alp[k] and wgt[k] (k>0) are meaningfull, 
       *     and it is illegal to reference alp[0] and wgt[0]
       */

      //propagation and analysis
      if(excited) {
        // It is in the excited state, propagation with 1/2 m w^2 (x-delx)^2
        //printf("        propagate at state = %4d\n", 0);
        p0   = fpt*imag(alp);
        alp  = (alp-dalp)*expfreqdt + dalp;
        pt   = fpt*imag(alp);
        wgt *= exphwdt*exp(0.5*I*(p0-pt)*delx);
        //generate the wavefunction for average
      } else {
        // It is in the continuum state, propagation with 1/2 m w^2 x^2
        //printf("        propagate at state = %4d\n", idx);
        wgt *= exphwdt;
        alp *= expfreqdt;
      }
      for(int l=0; l<Nel; l++) wgt *= expEndt[state[l]];
      jc[nj]++;

      if(j%nbloc==0) {
#ifdef _YYY_REALSPACE_AV_
        cswf(npt, xmin, dx, mass, freq, wgt*Iton[nj%4]*sclf[j], alp, wf[j/nbloc][istate]);
#else
        wfCS2En(Neig, wgt*Iton[nj%4]*sclf[j], alp, aux, wf[j/nbloc][istate]);
#endif
      }
    }
	}

#ifdef _YYY_REALSPACE_AV_
  int nmem = npt;
#else
  int nmem = Neig+1;
#endif

#ifdef _YYY_MPI_
  dcomplex ***avg  = array3d<dcomplex>(nwf+1, Nhs, nmem);
  MPI_Allreduce(**wf, **avg,  (1+nwf)*Nhs*nmem,  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  memcpy(**wf, **avg, sizeof(dcomplex)*(1+nwf)*Nhs*nmem);
  free3d(avg);

  int *jcl  = array1d<int>(nstep);
  MPI_Allreduce(jc, jcl,  nstep, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(jc, jcl, sizeof(int)*nstep);
  free1d(jcl);
#endif

  if(myid != master) {
    free2d(wf);
    free1d(sclf);
    free1d(jc);
    return;
  }

  char fnm[256];
  sprintf(fnm, "ahm-mbcs-s%d-n%d-%d.dat", Norb, Nel, ntraj);
  FILE *FL = fopen(fnm, "w");

  fprintf(FL, "#discretizing the bath\n");
  for(int n=0; n<Norb; n++) fprintf(FL, "#%6d %1.16e %1.16e\n", n, cpl[n], En[n]);
  for(int t=0; t<=nwf; t++)  wfanalysis(FL, t*dt*nbloc, wf[t]);
  fclose(FL);

  sprintf(fnm, "ahm-sepjc-s%d-n%d-%d.dat", Norb, Nel, ntraj);
  FL = fopen(fnm, "w");
  for(int t=0; t<nstep; t++) {
    fprintf(FL, "%20d %1.16e\n", t, jc[t]*1.0/ntraj);
  }
  fclose(FL);

  free2d(wf);
  free1d(sclf);
  free1d(jc);

	return;
}

