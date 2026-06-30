
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
  int nbloc = 1;
  double **prb = array2d<double>(nwf+1, Norb+3);
  int excited0 = S0.find(0) == S0.end() ? 0 : 1,
      Jmax     = nstep + 1;
  csproj(wgt_init*Lfct, alp_init, wgt_init, alp_init, excited0, S0, prb[0]);

  const double lambda = abs(cpl[1]),
               rate   = sqrtNel*sqrtNvac*lambda*dt;
  double p0, pt,  inv_ntraj = 1.0/ntraj,
         *sclf       = array1d<double>(nstep+1);
  int    *jumps_back = array1d<int>(Jmax),
         idx;
  for(int j=0; j<=nstep;  j++) sclf[j]       = exp(rate*j);

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
        if((int)state.size() != Nel || (int)vac.size() != Nvac ||
           (excited && state.find(0) == state.end()) ||
           (!excited && state.find(0) != state.end())) {
          cerr<<"invalid forward SepMB state before jump.\n";
          abort();
        }
        //idx in this block indicate the orbital involved in the quantum jump
        if(excited) {
          // if the molecule is excited, the electron jumps to an vacant orbital
          if(vac.empty()) {
            cerr<<"empty vacant set in excited SepMB jump.\n";
            abort();
          }
          idx = get_random_element(vac);
          if(idx < 0 || idx >= Norb) {
            cerr<<"invalid vacant orbital sampled in SepMB jump.\n";
            abort();
          }
          state.erase(0);
          state.insert(idx);
#ifdef _CHECK_PATH_
          path.push_back({0, idx});
#endif
          vac.insert(0);
          vac.erase(idx);
          int pos = 0;
          auto it = state.begin();
          for(; it != state.end() && *it != idx; ++it, ++pos) {}
          if(it == state.end()) {
            cerr<<"failed to locate excited forward SepMB jump target.\n";
            abort();
          }
          sign = eo[pos%2];
#ifdef _TRACE_STATE_
          printf("  forward %2d-th jump : sign = %+d, switch 0 -> %4d : \n", nj, sign, idx);
#endif
          wgt *= sign*sqrtNvac*sqrt1_Nel;
        } else {
          // if the molecule is not excited, the electron jumps from an occupied orbital to the molecule
          if(state.empty()) {
            cerr<<"empty occupied set in ground SepMB jump.\n";
            abort();
          }
          idx  = get_random_element(state);
          if(idx <= 0 || idx >= Norb) {
            cerr<<"invalid occupied orbital sampled in SepMB jump.\n";
            abort();
          }
          int pos = 0;
          auto it = state.begin();
          for(; it != state.end() && *it != idx; ++it, ++pos) {}
          if(it == state.end()) {
            cerr<<"failed to locate ground forward SepMB jump target.\n";
            abort();
          }
          sign = eo[pos%2];
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

      if(j<=nwf*nbloc && j%nbloc==0) {
        int excited_for = excited;
        // the distance between S0 and state in Johnson graph.
        set<int> C;
        set_difference(state.begin(), state.end(), S0.begin(), S0.end(), 
                   inserter(C, C.begin()));
        int d = C.size(),
            class_flip = excited0 ^ excited,
            d_kondo = d - class_flip,
            nj_min = 2*d_kondo + class_flip;
        if(d_kondo < 0 || d_kondo > min(Nel, Nvac)) {
          cerr<<"invalid Kondo distance for separated many-body projection.\n";
          abort();
        }
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
        set<int> state_for = state;
        int offset  = 0;
        alp     = alp_init;
        wgt     = wgt_init;
        excited = excited0;
        state   = S0;
        bool valid_path = ((int)path.size() == nj_back);
        // vac is no longer needed in backward propagation
        //set_difference(orbitals.begin(), orbitals.end(), state.begin(), state.end(), 
        //           inserter(vac, vac.begin()));
        for(int k=0; valid_path && k<nj_back; k++) {
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
          if(path[k].first < 0 || path[k].first >= Norb ||
             path[k].second < 0 || path[k].second >= Norb ||
             state.find(path[k].first) == state.end() ||
             state.find(path[k].second) != state.end()) {
            valid_path = false;
            break;
          }
          if(excited) {
            state.erase(path[k].first);
            state.insert(path[k].second);
            idx = path[k].second;
          } else {
            idx = path[k].first;
            state.erase(path[k].first);
            state.insert(path[k].second);
          }
          int pos = 0;
          auto it = state.begin();
          for(; it != state.end() && *it != idx; ++it, ++pos) {}
          if(it == state.end() || (int)state.size() != Nel) {
            valid_path = false;
            break;
          }
          sign = eo[pos%2];
#ifdef _TRACE_STATE_
          printf("  backward %2d-th jump for time %d, ", k, j);
          printf(   "     sign = %+d, switch %4d <-> %4d : \n", sign, path[k].first, path[k].second);
#endif
          wgt *= sign*sqrtfct[excited];
          excited = 1 - excited;
        }

        if(!valid_path) {
          alp = alp_for;
          wgt = wgt_for;
          state = state_for;
          excited = excited_for;
          continue;
        }

        // The odd Kondo sector is a virtual charge-transfer sector in the
        // separated estimator.  Counting it directly as a diagonal population
        // makes orbital 0 deplete and overpopulates high bath orbitals.
        // Match the exact many-body observable by projecting populations only
        // after the path has returned to the initial electronic class.
        if(excited_for != excited0) {
          alp = alp_for;
          wgt = wgt_for;
          state = state_for;
          excited = excited_for;
          continue;
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

        double measure = 1.0; // jumps_back is already sampled with a lambda-dependent Poisson rate
        double endpoint_prb = excited0 ? sampler.get_Ptd(nj_back, d_kondo)
                                       : sampler.get_Qtd(nj_back, d_kondo);
        if(endpoint_prb <= 0.0) {
          alp = alp_for;
          wgt = wgt_for;
          state = state_for;
          excited = excited_for;
          continue;
        }
        measure *= endpoint_prb*sclf[j-offset];
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
        state = state_for;
        excited = excited_for;
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
    fprintf(FL, "#PATCH_CHECK: SepMBpoisson v0.43 even-sector-observable-projection trace-normalized output active\n");
    fprintf(FL, "#discretizing the bath:\n");
    for(int n=0; n<Norb; n++) fprintf(FL, "#%6d %1.16e %1.16e\n", n, cpl[n], En[n]);
    for(int t=0; t<=nwf; t++)  {
      double norm = prb[t][0]/Nel;
      if(fabs(norm) > 1.e-300) {
        fprintf(FL, "%12.8f %+1.16e %+1.16e %+1.16e", t*dt, (double)Nel, prb[t][1]/norm, prb[t][2]/norm);
        for(int k=3; k<Norb+3; k++) fprintf(FL, " %+1.16e", prb[t][k]/norm);
      } else {
        fprintf(FL, "%12.8f %+1.16e %+1.16e %+1.16e", t*dt, prb[t][0], prb[t][1], prb[t][2]);
        for(int k=3; k<Norb+3; k++) fprintf(FL, " %+1.16e", prb[t][k]);
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

