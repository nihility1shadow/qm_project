
#include <gsl/gsl_sf_lambert.h>
#ifdef _YYY_MPI_
#include <mpi.h>
//#include "./yyy_mpi.h"
#endif

#include "./ahm.h"
#include "./qj.h"
#include "./qmfft.h"
#include "./qho.h"
#include "./Kondo-path-sampler.h"

#define _YYY_REALSPACE_AV_

/*
 * This file tries to obtain the average without direct construction
 * the final wavefunction.
 */

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
typedef int (KondoPathSampler::*Path_Sampler)(const int k, const set<int> &S0, const set<int> &S1,
              vector<pair<int, int> > &path) const;

  
void AHM::SepMBpoisson(const int ntraj, const int nstep, const double dt, 
    const dcomplex alp_init, const dcomplex wgt_init, const set<int> &S0) const {
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
              dalp     = mass*freq*delx/sqrt(2*mass*freq), // the shift of alpha parameter corresponding to delx
            sqrtfct[2] = {sqrtNel*sqrt1_Nvac , sqrtNvac*sqrt1_Nel};
  dcomplex *exphwdt    = array1d<dcomplex>(nstep+1), //exp(-0.5*freq*dt*I),
           *expfreqdt  = array1d<dcomplex>(nstep+1),
           **expEndt   = array2d<dcomplex>(nstep+1, Norb),
           *expEndt1, expfreqdt1, exphwdt1, alp, wgt;

  for(int j=2; j<Norb; j++) if(cpl[j-1] != cpl[j]) {
    printf("At this momentum, the Poisson jump method is only valid for homogeneous cj.\n");
    abort();
  }

  for(int k=0; k<=nstep; k++) {
    expfreqdt[k] = exp(-(k*freq*dt)*I);
    exphwdt[k]   = exp(-(0.5*freq*k*dt)*I);
    for(int j=0; j<Norb; j++) expEndt[k][j] = exp(-I*(dt*En[j]*k));
  }
  expEndt1   = expEndt[1];
  expfreqdt1 = expfreqdt[1];
  exphwdt1   = exphwdt[1];

  int ntraj_local = ntraj, *jc = array1d<int>(nstep+1);
#ifdef _YYY_MPI_
  ntraj_local = ntraj/nproc;
  if(myid==0) ntraj_local += ntraj - ntraj_local*nproc;
  const double Lfct = 1.0/nproc;
#else
  const double Lfct = 1.0;
#endif

  vector<pair<int, int>> path;
  set<int> state,    // the occupied orbitals
           orbitals, // all orbitals
           vac;      // the vacant orbitals
  for(int j=0; j<Norb; j++) orbitals.insert(j);

  // the initial average
  int nwf = nstep > 200 ? 200 : nstep;
  int nbloc = nstep/nwf;
  nwf = nstep/nbloc;
  double **prb = array2d<double>(nwf+1, Norb+3);
  int excited0 = S0.find(0) == S0.end() ? 0 : 1,
      Jmax     = nstep + 1;
  csproj(wgt_init*Lfct, alp_init, wgt_init, alp_init, excited0, S0, prb[0]);

  const double lambda = abs(cpl[1]),
               rate   = sqrtNel*sqrtNvac*lambda*dt;
  double p0, pt,  inv_ntraj = 1.0/ntraj,
         *pow_lambda = array1d<double>(Jmax),
         *sclf       = array1d<double>(nstep+1);
  int    *jumps_back = array1d<int>(Jmax),
         idx;
  pow_lambda[0] = 1.0;
  for(int j=0; j<=nstep;  j++) sclf[j]       = exp(rate*j);
  for(int j=1; j<Jmax; j++) pow_lambda[j] = pow_lambda[j-1]*lambda;

  KondoPathSampler sampler(Norb, Nel, Jmax);
  Path_Sampler PS[2][2] = {
          &KondoPathSampler::sample_path_B2B, 
          &KondoPathSampler::sample_path_B2A,
          &KondoPathSampler::sample_path_A2B,
          &KondoPathSampler::sample_path_A2A
  };

  for(int n=0; n<ntraj_local; n++) {
    alp    = alp_init;
    wgt    = wgt_init;
    state  = S0;
    vac.clear();
    set_difference(orbitals.begin(), orbitals.end(), state.begin(), state.end(), 
                   inserter(vac, vac.begin()));
#ifdef _TRACE_STATE_
        cout<<n<<"-th trajectory:\n";
        print_set("      occ", state, "\n");
        print_set("      vac", vac, "\n");
#endif
    int excited = excited0, //the index show the molecule is excited (1) or not (0)
        nj      = 0,        //number of jumps in the forward  path;
        sign, found;

#ifdef _CHECK_PATH_
    path.clear();
#endif
    //the forward propagation
    for(int j=1; j<=nstep; j++) {
      //printf("    time: %4d", j);
      // jump seperatedly, switch the state between 0 and idx
      if(drand48() < rate) {
        //idx in this block indicate the orbital involved in the quantum jump
        if(excited) {
          // if the molecule is excited, the electron jumps to an vacant orbital
          idx = get_random_element(vac);
          state.erase(0);
          state.insert(idx);
#ifdef _CHECK_PATH_
          path.push_back({0, idx});
#endif
          vac.insert(0);
          vac.erase(idx);
          auto it = state.begin();
          for(int l=0; l<Nel; l++,++it) if(*it == idx) {
            sign = eo[l%2];
            break;
          }
#ifdef _TRACE_STATE_
          printf("  forward %2d-th jump : sign = %+d, switch 0 -> %4d : \n", nj, sign, idx);
#endif
          wgt *= sign*sqrtNvac*sqrt1_Nel;
        } else {
          // if the molecule is not excited, the electron jumps from an occupied orbital to the molecule
          idx  = get_random_element(state);
          auto it = state.begin();
          for(int l=0; l<Nel; l++,++it) if(*it == idx) {
            sign = eo[l%2];
            break;
          }
          wgt *= sign*sqrtNel*sqrt1_Nvac;
#ifdef _TRACE_STATE_
          printf("  forward %2d-th jump : sign = %+d, switch 0 <- %4d : \n", nj, sign, idx);
#endif
          state.erase(idx);
          state.insert(0);
#ifdef _CHECK_PATH_
          path.push_back({idx, 0});
#endif
          vac.erase(0);
          vac.insert(idx);
        }
#ifdef _CHECK_PATH_
        jumps_back[nj] = j;
#endif
        nj++; 
        excited = 1 - excited;
#ifdef _TRACE_STATE_
//        cout<<"    "<<nj<<" jumps :\n";
//        print_set("      occ", state, "\n");
//        print_set("      vac", vac, "\n");
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
        alp  = (alp-dalp)*expfreqdt1 + dalp;
        pt   = fpt*imag(alp);
        wgt *= exphwdt1*exp(0.5*I*(p0-pt)*delx);
        //generate the wavefunction for average
      } else {
        // It is in the continuum state, propagation with 1/2 m w^2 x^2
        //printf("        propagate at state = %4d\n", idx);
        wgt *= exphwdt1;
        alp *= expfreqdt1;
      }
      for(int l : state) wgt *= expEndt1[l];
      jc[nj]++;

      if(j%nbloc==0) {
        int excited_for = excited;
        // the distance between S0 and state in Johnson graph.
        set<int> C;
        set_difference(state.begin(), state.end(), S0.begin(), S0.end(), 
                   inserter(C, C.begin()));
        int d = C.size(),
            nj_min = 2*d - (excited0 ^ excited);
#ifndef _CHECK_PATH_
        /* 
         * generate a Poisson process with nj_back jumps.
         * Note that: 
         *   1. nj_back should have the same parity as nj
         *   2. nj_back is sufficiently big to transform the system from S0 to state
         */
        int nj_back;
        do {
          nj_back = 0;
          bzero(jumps_back, Jmax*sizeof(int));
          for(int k=1; k<=j; k++)  if(drand48() < rate) jumps_back[nj_back++] = k;
        } while((nj_back+nj)%2 || nj_back<nj_min);

        int status = 0, fails=0;
        //print_set("state is: ", state, "\n");
        do {
          status = (sampler.*PS[excited0][excited])(nj_back, S0, state, path);
          fails -= status;
        } while(status<0 && fails < 5);
        if(status<0) continue;
#else
        int nj_back = nj;
        cout<<"    "<<nj<<" jumps in the forward path"<<endl;
        print_set("     state : ", state, "\n");
        cout<<"    the backward path, t = "<<j<<endl;
#endif

        //the backward paropagation, we should store the forward results first
        dcomplex alp_for = alp,
                 wgt_for = wgt;
        int offset  = 0;
        alp     = alp_init;
        wgt     = wgt_init;
        excited = excited0;
        state   = S0;
        // vac is no longer needed in backward propagation
        //set_difference(orbitals.begin(), orbitals.end(), state.begin(), state.end(), 
        //           inserter(vac, vac.begin()));
        for(int k=0; k<nj_back; k++) {
          //propagate the time offset -> jumps_back[k];
          int nadvance = jumps_back[k]-1-offset;
          if(excited) {
            // It is in the excited state, propagation with 1/2 m w^2 (x-delx)^2
            p0   = fpt*imag(alp);
            alp  = (alp-dalp)*expfreqdt[nadvance] + dalp;
            pt   = fpt*imag(alp);
            wgt *= exphwdt[nadvance]*exp(0.5*I*(p0-pt)*delx);
          } else {
            // It is in the continuum state, propagation with 1/2 m w^2 x^2
            wgt *= exphwdt[nadvance];
            alp *= expfreqdt[nadvance];
          }
          for(int l : state) wgt *= expEndt[nadvance][l];
          offset = jumps_back[k]-1;

          //deal with the k-th jump
          //idx in this block indicate the orbital involved in the quantum jump
          if(excited) {
            state.erase(path[k].first);
            state.insert(path[k].second);
            idx = path[k].second;
            auto it = state.begin();
            for(int l=0; l<Nel; l++,++it) if(*it == idx) {
              sign = eo[l%2];
              break;
            }
          } else {
            idx = path[k].first;
            auto it = state.begin();
            for(int l=0; l<Nel; l++,++it) if(*it == idx) {
              sign = eo[l%2];
              break;
            }
            state.erase(path[k].first);
            state.insert(path[k].second);
          }
#ifdef _TRACE_STATE_
          printf("  backward %2d-th jump for time %d, ", k, j);
          printf(   "     sign = %+d, switch %4d <-> %4d : \n", sign, path[k].first, path[k].second);
#endif
          wgt *= sign*sqrtfct[excited];
          excited = 1 - excited;
        }

        //propagate the time jump_bck[nj_back]->nstep;
        if(excited) {
          // It is in the excited state, propagation with 1/2 m w^2 (x-delx)^2
          //printf("        propagate at state = %4d\n", 0);
          p0   = fpt*imag(alp);
          alp  = (alp-dalp)*expfreqdt[j-offset] + dalp;
          pt   = fpt*imag(alp);
          wgt *= exphwdt[j-offset]*exp(0.5*I*(p0-pt)*delx);
        } else {
          // It is in the continuum state, propagation with 1/2 m w^2 x^2
          //printf("        propagate at state = %4d\n", idx);
          wgt *= exphwdt[j-offset];
          alp *= expfreqdt[j-offset];
        }
        for(int l : state) wgt *= expEndt[j-offset][l];

        if(!offset) offset = j;
        double measure = pow_lambda[nj_back]; // the probability density is lambda^nj*exp(-lambda*t)
        measure *= sampler.get_Ptd(nj_back, d)*sclf[j-offset];
        if(excited_for != excited) {
          cerr<<"parity mismatch between the forward and backward path.\n";
          abort();
        }
#ifdef _CHECK_PATH_
        print_set("     state : ", state, "\n");
        printf("wgt : %+1.16e %+1.16e %+1.16e %+1.16e\n", real(wgt_for), imag(wgt_for), real(wgt), imag(wgt));
        printf("alp : %+1.16e %+1.16e %+1.16e %+1.16e\n", real(alp_for), imag(alp_for), real(alp), imag(alp));
#endif
        csproj(wgt_for*Iton[nj%4]*sclf[j]*inv_ntraj, alp_for, wgt*Iton[nj_back%4]*measure, alp, excited, state, prb[j/nbloc]);

        //restore the forward parameters
        alp = alp_for;
        wgt = wgt_for;
      }
    }
	}

#ifdef _YYY_MPI_
  double **avg  = array2d<double>(nwf+1, Norb+3);
  MPI_Allreduce(*prb, *avg,  (1+nwf)*(Norb+3),  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  memcpy(*prb, *avg, sizeof(double)*(1+nwf)*(Norb+3));
  free2d(avg);

  int *jcl  = array1d<int>(nstep+1);
  MPI_Allreduce(jc, jcl,  nstep+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(jc, jcl, sizeof(int)*(nstep+1));
  free1d(jcl);
#endif

  if(myid==master) {
    char fnm[256];
    sprintf(fnm, "ahm-sepmb-s%d-n%d-%d.dat", Norb, Nel, ntraj);
    FILE *FL = fopen(fnm, "w");
    fprintf(FL, "#PATCH_CHECK: SepMBpoisson v0.32 fast-dense interpolated trace-normalized output active\n");
    fprintf(FL, "#discretizing the bath:\n");
    for(int n=0; n<Norb; n++) fprintf(FL, "#%6d %1.16e %1.16e\n", n, cpl[n], En[n]);
    for(int j=0; j<=nstep; j++)  {
      int t0 = j/nbloc;
      if(t0 > nwf) t0 = nwf;
      int t1 = t0 < nwf ? t0 + 1 : t0;
      double frac  = t1 == t0 ? 0.0 : (j - t0*nbloc)/(double)nbloc,
             norm0 = prb[t0][0]/Nel,
             norm1 = prb[t1][0]/Nel,
             v0, v1;
      v0 = fabs(norm0) > 1.e-300 ? (double)Nel : prb[t0][0];
      v1 = fabs(norm1) > 1.e-300 ? (double)Nel : prb[t1][0];
      fprintf(FL, "%12.8f %+1.16e", j*dt, (1.0-frac)*v0 + frac*v1);
      for(int k=1; k<Norb+3; k++) {
        v0 = fabs(norm0) > 1.e-300 ? prb[t0][k]/norm0 : prb[t0][k];
        v1 = fabs(norm1) > 1.e-300 ? prb[t1][k]/norm1 : prb[t1][k];
        fprintf(FL, " %+1.16e", (1.0-frac)*v0 + frac*v1);
      }
      fprintf(FL, "\n");
    }
    fclose(FL);

    sprintf(fnm, "ahm-jcmb-s%d-n%d-%d.dat", Norb, Nel, ntraj);
    FL = fopen(fnm, "w");
    for(int t=0; t<nstep; t++) {
      fprintf(FL, "%20d %1.16e\n", t, jc[t]*1.0/ntraj);
    }
    fclose(FL);
  }

  free2d(prb);
  free1d(sclf);
  free1d(jc);

	return;
}

