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
#include <cmath>
#include <cstdlib>
#include <map>

#define _YYY_REALSPACE_AV_

namespace {

double sepmb_binom(const int n, const int k) {
  if(k < 0 || k > n) return 0.0;
  int kk = k < n-k ? k : n-k;
  double r = 1.0;
  for(int j=1; j<=kk; j++) r *= (n-kk+j)*1.0/j;
  return r;
}

double sepmb_kondo_degeneracy(const int Norb, const int Nel,
    const int excited0, const int nj, const int d) {
  const int Nvac = Norb - Nel;
  if(excited0) {
    return (nj%2 == 0)
      ? sepmb_binom(Nel-1, d)*sepmb_binom(Nvac, d)
      : sepmb_binom(Nel-1, d)*sepmb_binom(Nvac, d+1);
  }
  return (nj%2 == 0)
    ? sepmb_binom(Nel, d)*sepmb_binom(Nvac-1, d)
    : sepmb_binom(Nel, d+1)*sepmb_binom(Nvac-1, d);
}

int sepmb_sample_conditioned_count(const double *tail, const int nstep,
    const int parity, const int kmin, const double u) {
  int first = kmin;
  if((first&1) != parity) first++;
  if(first > nstep || tail[first] <= 0.0) return -1;

  double target = u*tail[first];
  for(int k=first; k<=nstep; k+=2) {
    const double next_tail = k+2 <= nstep ? tail[k+2] : 0.0;
    const double mass = tail[k] - next_tail;
    if(target < mass || k+2 > nstep) return k;
    target -= mass;
  }
  return -1;
}

int sepmb_sample_jump_times(const int nstep, const int njump, int *jumps) {
  int selected = 0, remaining = njump;
  for(int j=1; j<=nstep && remaining>0; j++) {
    const int slots = nstep-j+1;
    if(drand48()*slots < remaining) {
      jumps[selected++] = j;
      remaining--;
    }
  }
  return remaining == 0 ? selected : -1;
}

int sepmb_sample_jump_times_stratified(const int nstep, const int njump,
    const unsigned long long sample_index, const unsigned long long shift,
    const bool stratify_single_jump, const bool stratify_multi_jump,
    int *jumps) {
  if(stratify_single_jump && njump == 1 && nstep > 0) {
    jumps[0] = 1 + (int)((sample_index+shift)%(unsigned long long)nstep);
    return 1;
  }
  if(stratify_multi_jump && njump > 1 && njump <= nstep) {
    // A random shift makes every sample unbiased.  The irrational rotation
    // spreads successive samples over the uniformly weighted combination
    // intervals, while this inverse map preserves the exact law for choosing
    // njump distinct grid points from nstep.
    const double inv_2_53 = 1.0/9007199254740992.0,
                 shift_u = (double)(shift>>11)*inv_2_53,
                 golden = 0.6180339887498948482;
    double u = fmod(shift_u + (sample_index+0.5)*golden, 1.0);
    int selected = 0, remaining = njump;
    for(int j=1; j<=nstep && remaining>0; j++) {
      const int slots = nstep-j+1;
      if(remaining == slots) {
        while(j<=nstep) jumps[selected++] = j++;
        remaining = 0;
        break;
      }
      const double include_probability = remaining/(double)slots;
      if(u < include_probability) {
        jumps[selected++] = j;
        remaining--;
        u /= include_probability;
      } else {
        u = (u-include_probability)/(1.0-include_probability);
      }
      if(u < 0.0) u = 0.0;
      if(u >= 1.0) u = nextafter(1.0, 0.0);
    }
    return remaining == 0 ? selected : -1;
  }
  return sepmb_sample_jump_times(nstep, njump, jumps);
}

unsigned long long sepmb_random_shift() {
  return ((unsigned long long)lrand48()<<32) ^ (unsigned long long)lrand48();
}

bool sepmb_resample_window_ensemble(const int nparticle, const int ndim,
    vector<dcomplex> &electronic, vector<dcomplex> &alpha,
    vector<dcomplex> &weight, vector<int> &excited,
    const double alpha_bin_width) {
  if(nparticle <= 0 || ndim <= 0) return false;

  vector<dcomplex> grouped_electronic;
  vector<dcomplex> grouped_alpha_sum;
  vector<double> grouped_alpha_weight;
  vector<int> grouped_excited;
  map<pair<int, pair<long long, long long> >, int> group_lookup;
  for(int n=0; n<nparticle; n++) {
    double norm2 = 0.0;
    const dcomplex *vec = &electronic[(size_t)n*ndim];
    for(int s=0; s<ndim; s++) norm2 += std::norm(vec[s]);
    const double particle_importance = abs(weight[n])*sqrt(norm2);
    if(!(particle_importance > 1.e-300)) continue;

    int group = -1;
    if(alpha_bin_width > 0.0) {
      const long long real_bin = (long long)floor(real(alpha[n])/alpha_bin_width),
                      imag_bin = (long long)floor(imag(alpha[n])/alpha_bin_width);
      const pair<int, pair<long long, long long> > key(
          excited[n], make_pair(real_bin, imag_bin));
      map<pair<int, pair<long long, long long> >, int>::iterator found =
          group_lookup.find(key);
      if(found == group_lookup.end()) {
        group = grouped_excited.size();
        group_lookup[key] = group;
      } else {
        group = found->second;
      }
    } else {
      group = grouped_excited.size();
    }

    if(group == (int)grouped_excited.size()) {
      grouped_excited.push_back(excited[n]);
      grouped_alpha_sum.push_back(0.0);
      grouped_alpha_weight.push_back(0.0);
      grouped_electronic.resize((size_t)(group+1)*ndim, 0.0);
    }
    dcomplex *group_vec = &grouped_electronic[(size_t)group*ndim];
    for(int s=0; s<ndim; s++) group_vec[s] += weight[n]*vec[s];
    grouped_alpha_sum[group] += particle_importance*alpha[n];
    grouped_alpha_weight[group] += particle_importance;
  }

  const int ngroup = grouped_excited.size();
  if(ngroup <= 0) return false;
  vector<double> group_norm(ngroup, 0.0), cumulative(ngroup, 0.0);
  double total = 0.0;
  for(int group=0; group<ngroup; group++) {
    double norm2 = 0.0;
    const dcomplex *vec = &grouped_electronic[(size_t)group*ndim];
    for(int s=0; s<ndim; s++) norm2 += std::norm(vec[s]);
    group_norm[group] = sqrt(norm2);
    total += group_norm[group];
    cumulative[group] = total;
  }
  if(!(total > 1.e-300) || !std::isfinite(total)) return false;

  vector<dcomplex> electronic_new((size_t)nparticle*ndim),
                   alpha_new(nparticle), weight_new(nparticle);
  vector<int> excited_new(nparticle);
  const double common_magnitude = total/nparticle,
               u0 = drand48()/nparticle;
  int source = 0;
  for(int n=0; n<nparticle; n++) {
    const double target = (u0+n/(double)nparticle)*total;
    while(source+1<ngroup && cumulative[source] < target) source++;
    if(group_norm[source] <= 1.e-300 ||
       grouped_alpha_weight[source] <= 1.e-300) return false;
    weight_new[n] = common_magnitude;
    alpha_new[n] = grouped_alpha_sum[source]/grouped_alpha_weight[source];
    excited_new[n] = grouped_excited[source];
    const dcomplex *src = &grouped_electronic[(size_t)source*ndim];
    dcomplex *dst = &electronic_new[(size_t)n*ndim];
    for(int s=0; s<ndim; s++) dst[s] = src[s]/group_norm[source];
  }

  electronic.swap(electronic_new);
  alpha.swap(alpha_new);
  weight.swap(weight_new);
  excited.swap(excited_new);
  return true;
}

vector<double> sepmb_binomial_cdf(const int nstep, const double probability) {
  vector<double> cdf(nstep+1, 0.0);
  long double mass = powl(1.0-probability, nstep);
  long double total = mass;
  cdf[0] = (double)total;
  for(int k=1; k<=nstep; k++) {
    mass *= (nstep-k+1)*probability/(k*(1.0-probability));
    total += mass;
    cdf[k] = (double)total;
  }
  if(total <= 0.0) return cdf;
  for(double &value : cdf) value /= (double)total;
  cdf.back() = 1.0;
  return cdf;
}

int sepmb_sample_cdf(const vector<double> &cdf, const double u) {
  const auto pos = lower_bound(cdf.begin(), cdf.end(), u);
  return pos == cdf.end() ? (int)cdf.size()-1 : (int)(pos-cdf.begin());
}


}

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
  if(myid==master && Nhs == Norb) {
    const int oldqm_nwf = nstep > 200 ? 200 : nstep;
    const double oldqm_lambda = abs(cpl[1]),
                 oldqm_sqrtN = sqrt((double)(Norb-1)),
                 oldqm_sqrt1_N = 1.0/oldqm_sqrtN,
                 oldqm_rate = oldqm_sqrtN*oldqm_lambda*dt;
    dcomplex **amp = array2d<dcomplex>(oldqm_nwf+1, Norb),
             *expE = array1d<dcomplex>(Norb);
    double *sclf_oldqm = array1d<double>(nstep+1);
    for(int j=0; j<=nstep; j++) sclf_oldqm[j] = exp(oldqm_rate*j)/ntraj;
    for(int n=0; n<Norb; n++) expE[n] = exp(-I*dt*Eocc[n]);
    amp[0][0] = wgt_init;

    vector<set<int> > qm_state(Norb);
    for(int n=0; n<Norb; n++) {
      for(int k=0; k<Nel; k++) qm_state[n].insert(occ[n][k]);
    }

    int *jc_oldqm = array1d<int>(nstep+1);
    for(int n=0; n<ntraj; n++) {
      int idx = 0, nj = 0;
      dcomplex wgt = wgt_init;
      for(int j=1; j<=nstep; j++) {
        if(drand48() < oldqm_rate) {
          nj++;
          if(idx) {
            wgt *= oldqm_sqrt1_N;
            idx = 0;
          } else {
            idx = (int)(drand48()*(Norb-1)) + 1;
            if(idx >= Norb) idx = Norb-1;
            wgt *= oldqm_sqrtN;
          }
        }
        wgt *= expE[idx];
        if(j <= oldqm_nwf) amp[j][idx] += wgt*Iton[nj%4]*sclf_oldqm[j];
        if(nj <= nstep) jc_oldqm[nj]++;
      }
    }

    char fnm[256];
    sprintf(fnm, "ahm-sepmb-s%d-n%d-%d.dat", Norb, Nel, ntraj);
    FILE *FL = fopen(fnm, "w");
    fprintf(FL, "#PATCH_CHECK: SepMBpoisson v0.52 old-qm-star-poisson active\n");
    fprintf(FL, "#discretizing the bath:\n");
    for(int n=0; n<Norb; n++) fprintf(FL, "#%6d %1.16e %1.16e\n", n, cpl[n], En[n]);
    for(int t=0; t<=oldqm_nwf; t++) {
      double *rlt = array1d<double>(Norb+3);
      for(int n=0; n<Norb; n++) {
        csproj(amp[t][n], alp_init, amp[t][n], alp_init, 1, qm_state[n], rlt);
      }
      double norm = rlt[0]/Nel;
      if(fabs(norm) > 1.e-300) {
        fprintf(FL, "%12.8f %+1.16e %+1.16e %+1.16e", t*dt, (double)Nel, rlt[1]/norm, rlt[2]/norm);
        for(int k=3; k<Norb+3; k++) fprintf(FL, " %+1.16e", rlt[k]/norm);
      } else {
        fprintf(FL, "%12.8f %+1.16e %+1.16e %+1.16e", t*dt, rlt[0], rlt[1], rlt[2]);
        for(int k=3; k<Norb+3; k++) fprintf(FL, " %+1.16e", rlt[k]);
      }
      fprintf(FL, "\n");
      free1d(rlt);
    }
    fclose(FL);

    sprintf(fnm, "ahm-jcmb-s%d-n%d-%d.dat", Norb, Nel, ntraj);
    FL = fopen(fnm, "w");
    for(int t=0; t<nstep; t++) fprintf(FL, "%20d %1.16e\n", t, jc_oldqm[t]*1.0/ntraj);
    fclose(FL);

    free2d(amp);
    free1d(expE);
    free1d(sclf_oldqm);
    free1d(jc_oldqm);
    return;
  }

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

  int ntraj_local = ntraj, trajectory_offset = 0,
      *jc = array1d<int>(nstep+1);
#ifdef _YYY_MPI_
  const int trajectories_per_rank = ntraj/nproc,
            trajectory_remainder = ntraj-trajectories_per_rank*nproc;
  ntraj_local = trajectories_per_rank;
  if(myid==0) ntraj_local += trajectory_remainder;
  trajectory_offset = myid == 0 ? 0
      : trajectory_remainder + myid*trajectories_per_rank;
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
  const char *env_tmax = getenv("SEP_MB_TMAX"),
             *env_nwf  = getenv("SEP_MB_NWF"),
             *env_measure_stride = getenv("SEP_MB_MEASURE_STRIDE"),
             *env_back_replicas = getenv("SEP_MB_BACK_REPLICAS"),
             *env_stratify_forward = getenv("SEP_MB_STRATIFY_FORWARD_COUNT"),
             *env_exact_orbitals = getenv("SEP_MB_EXACT_ORBITALS"),
             *env_stratify_single_jump = getenv("SEP_MB_STRATIFY_SINGLE_JUMP_TIME"),
             *env_stratify_multi_jump = getenv("SEP_MB_STRATIFY_MULTI_JUMP_TIME"),
             *env_sample_back_orbitals = getenv("SEP_MB_SAMPLE_BACK_ORBITALS"),
             *env_exact_back_jumps = getenv("SEP_MB_EXACT_BACK_JUMPS"),
             *env_window_steps = getenv("SEP_MB_WINDOW_STEPS"),
             *env_window_pair_replicas = getenv("SEP_MB_WINDOW_PAIR_REPLICAS"),
             *env_window_alpha_bin = getenv("SEP_MB_WINDOW_ALPHA_BIN");
  double output_tmax = env_tmax ? atof(env_tmax) : 500.0;
  int nwf = output_tmax > 0.0 ? (int)(output_tmax/dt + 0.5) : 1000;
  if(env_nwf) nwf = atoi(env_nwf);
  if(nwf < 0) nwf = 0;
  if(nwf > nstep) nwf = nstep;
  int forced_measure_stride = env_measure_stride ? atoi(env_measure_stride) : 1;
  if(forced_measure_stride < 0) forced_measure_stride = 0;
  int back_replicas = env_back_replicas ? atoi(env_back_replicas) : 256;
  if(back_replicas < 1) back_replicas = 1;
  if(back_replicas > 1024) back_replicas = 1024;
  const bool stratify_forward = env_stratify_forward
      ? atoi(env_stratify_forward) != 0 : true;
  const bool exact_orbitals = env_exact_orbitals
      ? atoi(env_exact_orbitals) != 0 : true;
  const bool stratify_single_jump = env_stratify_single_jump
      ? atoi(env_stratify_single_jump) != 0 : true;
  const bool stratify_multi_jump = env_stratify_multi_jump
      ? atoi(env_stratify_multi_jump) != 0 : false;
  const bool sample_back_orbitals = env_sample_back_orbitals
      ? atoi(env_sample_back_orbitals) != 0 : true;
  int exact_back_jumps = env_exact_back_jumps ? atoi(env_exact_back_jumps) : 1;
  if(exact_back_jumps < 0) exact_back_jumps = 0;
  if(exact_back_jumps > 6) exact_back_jumps = 6;
  int window_steps = env_window_steps ? atoi(env_window_steps) : 0;
  if(window_steps < 0) window_steps = 0;
  if(window_steps > nwf) window_steps = nwf;
  int window_pair_replicas = env_window_pair_replicas
      ? atoi(env_window_pair_replicas) : 1;
  if(window_pair_replicas < 1) window_pair_replicas = 1;
  if(window_pair_replicas > 64) window_pair_replicas = 64;
  double window_alpha_bin = env_window_alpha_bin
      ? atof(env_window_alpha_bin) : 0.0;
  if(window_alpha_bin < 0.0) window_alpha_bin = 0.0;
  double gap = 0.0;
  for(int a : S0) {
    for(int b=0; b<Norb; b++) {
      if(S0.find(b) != S0.end()) continue;
      double de = fabs(En[a]-En[b]);
      if(de > 1.e-12 && (gap == 0.0 || de < gap)) gap = de;
    }
  }
  if(gap == 0.0) {
    for(int a=0; a<Norb; a++) {
      for(int b=a+1; b<Norb; b++) {
        double de = fabs(En[a]-En[b]);
        if(de > 1.e-12 && (gap == 0.0 || de < gap)) gap = de;
      }
    }
  }
  const int period_steps = gap > 0.0 ? max(1, (int)(2.0*acos(-1.0)/(gap*dt)+0.999999999999)) : nwf;
  const int dense_end    = min(nwf, max(16, period_steps/8));
  const int mid_end      = min(nwf, max(dense_end, period_steps/4));
  const int slow_end     = min(nwf, max(mid_end, period_steps/2));
  const int stride_mid   = max(1, period_steps/128);
  const int stride_slow  = max(stride_mid, period_steps/64);
  const int stride_late  = max(stride_slow, period_steps/32);
  vector<int> measure_steps;
  measure_steps.push_back(0);
  int last_step = 0;
  for(int j=1; j<=nwf; ) {
    int stride = forced_measure_stride > 0 ? forced_measure_stride :
                 (j <= dense_end ? 1 : (j <= mid_end ? stride_mid :
                 (j <= slow_end ? stride_slow : stride_late)));
    if(j > last_step) {
      measure_steps.push_back(j);
      last_step = j;
    }
    j += stride;
  }
  if(measure_steps.back() != nwf) measure_steps.push_back(nwf);
  int nmeas = measure_steps.size()-1;
  int *measure_slot = array1d<int>(nwf+1);
  for(int j=0; j<=nwf; j++) measure_slot[j] = -1;
  for(int j=0; j<=nmeas; j++) measure_slot[measure_steps[j]] = j;
  double **prb = array2d<double>(nmeas+1, Norb+3);
  int excited0 = S0.find(0) == S0.end() ? 0 : 1,
      Jmax     = nstep + 1;
  csproj(wgt_init*Lfct, alp_init, wgt_init, alp_init, excited0, S0, prb[0]);

  const double lambda          = abs(cpl[1]),
               jump_strength   = sqrtNel*sqrtNvac*lambda*dt,
               jump_probability = jump_strength/(1.0+jump_strength),
               log_scale       = log1p(jump_strength),
               inv_back_replicas = 1.0/back_replicas,
               inv_jump_normalization = 1.0/(lambda*sqrtNel*sqrtNvac);
  double p0, pt,  inv_ntraj = 1.0/ntraj,
          *sclf       = array1d<double>(nstep+1);
  int    *jumps_back = array1d<int>(Jmax),
         *jumps_forward = array1d<int>(Jmax),
         *forward_jump_schedule = array1d<int>(Jmax),
         idx;
  for(int j=0; j<=nstep;  j++) sclf[j] = exp(log_scale*j);

  const vector<double> forward_count_cdf = stratify_forward
      ? sepmb_binomial_cdf(nstep, jump_probability) : vector<double>();
  double forward_count_shift = 0.0;
  if(stratify_forward) {
#ifdef _YYY_MPI_
    if(myid==master) forward_count_shift = drand48();
    MPI_Bcast(&forward_count_shift, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
#else
    forward_count_shift = drand48();
#endif
  }
  unsigned long long forward_time_shift = 0,
                     backward_time_shift = 0;
  if(stratify_single_jump || stratify_multi_jump) {
#ifdef _YYY_MPI_
    if(myid==master) {
      forward_time_shift = sepmb_random_shift();
      backward_time_shift = sepmb_random_shift();
    }
    MPI_Bcast(&forward_time_shift, 1, MPI_UNSIGNED_LONG_LONG, master, MPI_COMM_WORLD);
    MPI_Bcast(&backward_time_shift, 1, MPI_UNSIGNED_LONG_LONG, master, MPI_COMM_WORLD);
#else
    forward_time_shift = sepmb_random_shift();
    backward_time_shift = sepmb_random_shift();
#endif
  }

  double ***back_accept = array3d<double>(nwf+1, 2, nwf+2);
  for(int jt=0; jt<=nwf; jt++) {
    double *pk = array1d<double>(jt+1);
    pk[0] = pow(1.0-jump_probability, jt);
    for(int k=1; k<=jt; k++) {
      pk[k] = pk[k-1]*(jt-k+1)*jump_probability/
              (k*(1.0-jump_probability));
    }
    for(int parity=0; parity<2; parity++) {
      double tail = 0.0;
      back_accept[jt][parity][jt+1] = 0.0;
      for(int k=jt; k>=0; k--) {
        if((k&1) == parity) tail += pk[k];
        back_accept[jt][parity][k] = tail;
      }
    }
    free1d(pk);
  }

  if(exact_orbitals) {
    if(lambda <= 0.0) {
      cerr<<"exact-orbital SepMB requires a nonzero homogeneous coupling.\n";
      abort();
    }

    vector<set<int> > basis_states(Nhs);
    int initial_basis = -1;
    for(int s=0; s<Nhs; s++) {
      for(int k=0; k<Nel; k++) basis_states[s].insert(occ[s][k]);
      if(basis_states[s] == S0) initial_basis = s;
    }
    if(initial_basis < 0) {
      cerr<<"failed to locate the initial determinant in SepMB basis.\n";
      abort();
    }

    dcomplex **expEoccdt = array2d<dcomplex>(nstep+1, Nhs),
             *vec_for    = array1d<dcomplex>(Nhs),
             *vec_back   = sample_back_orbitals ? NULL : array1d<dcomplex>(Nhs),
             *vec_tmp    = array1d<dcomplex>(Nhs);
    const int vec_bytes = sizeof(dcomplex)*Nhs;
    for(int k=0; k<=nstep; k++) {
      for(int s=0; s<Nhs; s++) {
        expEoccdt[k][s] = exp(-I*(dt*Eocc[s]*(double)k));
      }
    }

    if(window_steps > 0) {
      if(ntraj_local <= 0) {
        cerr<<"windowed SepMB requires at least one trajectory per MPI rank.\n";
        abort();
      }

      const size_t ensemble_size = (size_t)ntraj_local*Nhs;
      vector<dcomplex> ket_electronic(ensemble_size, 0.0),
                       bra_electronic(ensemble_size, 0.0),
                       ket_alpha(ntraj_local, alp_init),
                       bra_alpha(ntraj_local, alp_init),
                       ket_weight(ntraj_local, wgt_init),
                       bra_weight(ntraj_local, wgt_init);
      vector<int> ket_excited(ntraj_local, excited0),
                  bra_excited(ntraj_local, excited0),
                  ket_njump(ntraj_local, 0),
                  bra_njump(ntraj_local, 0);
      vector<unsigned char> ket_schedule((size_t)ntraj_local*(window_steps+1), 0),
                            bra_schedule((size_t)ntraj_local*(window_steps+1), 0);
      for(int n=0; n<ntraj_local; n++) {
        ket_electronic[(size_t)n*Nhs+initial_basis] = 1.0;
        bra_electronic[(size_t)n*Nhs+initial_basis] = 1.0;
      }

      const auto prepare_window_schedule = [&](vector<unsigned char> &schedule,
          const int block_steps, const int block_index,
          const double count_shift, const unsigned long long time_shift,
          int *jump_buffer) {
        fill(schedule.begin(), schedule.end(), 0);
        const vector<double> count_cdf = sepmb_binomial_cdf(
            block_steps, jump_probability);
        for(int n=0; n<ntraj_local; n++) {
          const unsigned long long global_index =
              (unsigned long long)block_index*ntraj + trajectory_offset+n;
          const double count_u = fmod(count_shift +
              (global_index+0.5)*0.6180339887498948482, 1.0);
          const int count = sepmb_sample_cdf(count_cdf, count_u);
          if(sepmb_sample_jump_times_stratified(block_steps, count,
                global_index, time_shift, stratify_single_jump,
                stratify_multi_jump, jump_buffer) != count) {
            cerr<<"failed to sample a windowed SepMB jump schedule.\n";
            abort();
          }
          for(int k=0; k<count; k++) {
            schedule[(size_t)n*(window_steps+1)+jump_buffer[k]] = 1;
          }
        }
      };

      const auto advance_window_ensemble = [&](vector<dcomplex> &electronic,
          vector<dcomplex> &alpha, vector<dcomplex> &weight,
          vector<int> &excited, vector<int> &njump,
          const vector<unsigned char> &schedule, const int local_step) {
        for(int n=0; n<ntraj_local; n++) {
          dcomplex *vec = &electronic[(size_t)n*Nhs];
          if(schedule[(size_t)n*(window_steps+1)+local_step]) {
            bzero(vec_tmp, vec_bytes);
            exc.multiply(vec, vec_tmp);
            for(int s=0; s<Nhs; s++) {
              vec[s] = vec_tmp[s]*inv_jump_normalization;
            }
            njump[n]++;
            excited[n] = 1-excited[n];
          }

          if(excited[n]) {
            const double p_before = fpt*imag(alpha[n]);
            alpha[n] = (alpha[n]-dalp)*expfreqdt1 + dalp;
            const double p_after = fpt*imag(alpha[n]);
            weight[n] *= exphwdt1
                *exp(0.5*I*(p_before-p_after)*delx);
          } else {
            weight[n] *= exphwdt1;
            alpha[n] *= expfreqdt1;
          }
          for(int s=0; s<Nhs; s++) vec[s] *= expEoccdt[1][s];
        }
      };

      int block_index = 0;
      for(int block_start=0; block_start<nwf;
          block_start += window_steps, block_index++) {
        const int block_steps = min(window_steps, nwf-block_start);
        double ket_count_shift = 0.0, bra_count_shift = 0.0;
        unsigned long long ket_time_shift = 0, bra_time_shift = 0;
#ifdef _YYY_MPI_
        if(myid==master) {
          ket_count_shift = drand48();
          bra_count_shift = drand48();
          ket_time_shift = sepmb_random_shift();
          bra_time_shift = sepmb_random_shift();
        }
        MPI_Bcast(&ket_count_shift, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&bra_count_shift, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
        MPI_Bcast(&ket_time_shift, 1, MPI_UNSIGNED_LONG_LONG, master, MPI_COMM_WORLD);
        MPI_Bcast(&bra_time_shift, 1, MPI_UNSIGNED_LONG_LONG, master, MPI_COMM_WORLD);
#else
        ket_count_shift = drand48();
        bra_count_shift = drand48();
        ket_time_shift = sepmb_random_shift();
        bra_time_shift = sepmb_random_shift();
#endif
        prepare_window_schedule(ket_schedule, block_steps, block_index,
            ket_count_shift, ket_time_shift, jumps_forward);
        prepare_window_schedule(bra_schedule, block_steps, block_index,
            bra_count_shift, bra_time_shift, jumps_back);
        fill(ket_njump.begin(), ket_njump.end(), 0);
        fill(bra_njump.begin(), bra_njump.end(), 0);

        for(int local_step=1; local_step<=block_steps; local_step++) {
          advance_window_ensemble(ket_electronic, ket_alpha, ket_weight,
              ket_excited, ket_njump, ket_schedule, local_step);
          advance_window_ensemble(bra_electronic, bra_alpha, bra_weight,
              bra_excited, bra_njump, bra_schedule, local_step);

          const int global_step = block_start+local_step,
                    iprb = measure_slot[global_step];
          for(int n=0; n<ntraj_local; n++) {
            if(ket_njump[n] <= nstep) jc[ket_njump[n]]++;
          }
          if(iprb < 0) continue;

          vector<int> bra_class[2];
          bra_class[0].reserve(ntraj_local);
          bra_class[1].reserve(ntraj_local);
          for(int m=0; m<ntraj_local; m++) {
            bra_class[bra_excited[m]&1].push_back(m);
          }
          int pair_offset[2] = {0, 0};
          for(int c=0; c<2; c++) {
            if(!bra_class[c].empty()) {
              pair_offset[c] = (int)(drand48()*bra_class[c].size());
              if(pair_offset[c] >= (int)bra_class[c].size()) {
                pair_offset[c] = bra_class[c].size()-1;
              }
            }
          }
          for(int n=0; n<ntraj_local; n++) {
            const int c = ket_excited[n]&1,
                      class_size = bra_class[c].size();
            if(class_size <= 0) continue;
            const double class_fraction = class_size/(double)ntraj_local,
                         inv_pair_replicas = 1.0/window_pair_replicas;
            const dcomplex ket_scale = ket_weight[n]
                *Iton[ket_njump[n]%4]*sclf[local_step]*inv_ntraj
                *class_fraction*inv_pair_replicas;
            const dcomplex *ket_vec = &ket_electronic[(size_t)n*Nhs];
            const int pair_stride = max(1, class_size/window_pair_replicas);
            for(int replica=0; replica<window_pair_replicas; replica++) {
              const int pair_position = (pair_offset[c]+n+
                      replica*pair_stride)%class_size,
                        m = bra_class[c][pair_position];
              const dcomplex bra_scale = bra_weight[m]
                  *Iton[bra_njump[m]%4]*sclf[local_step];
              const dcomplex *bra_vec = &bra_electronic[(size_t)m*Nhs];
              for(int s=0; s<Nhs; s++) {
                if(abs(ket_vec[s]) <= 1.e-300 ||
                   abs(bra_vec[s]) <= 1.e-300) continue;
                csproj(ket_scale*ket_vec[s], ket_alpha[n],
                       bra_scale*bra_vec[s], bra_alpha[m],
                       ket_excited[n], basis_states[s], prb[iprb]);
              }
            }
          }
        }

        for(int n=0; n<ntraj_local; n++) {
          ket_weight[n] *= Iton[ket_njump[n]%4]*sclf[block_steps];
          bra_weight[n] *= Iton[bra_njump[n]%4]*sclf[block_steps];
        }
        if(block_start+block_steps < nwf) {
          if(!sepmb_resample_window_ensemble(ntraj_local, Nhs,
                ket_electronic, ket_alpha, ket_weight, ket_excited,
                window_alpha_bin) ||
             !sepmb_resample_window_ensemble(ntraj_local, Nhs,
                bra_electronic, bra_alpha, bra_weight, bra_excited,
                window_alpha_bin)) {
            cerr<<"windowed SepMB resampling failed.\n";
            abort();
          }
        }
      }
      free2d(expEoccdt);
      free1d(vec_for);
      if(vec_back) free1d(vec_back);
      free1d(vec_tmp);
    } else {
    vector<vector<int> > back_orbital_targets;
    vector<vector<double> > back_orbital_factors;
    dcomplex *sparse_back = sample_back_orbitals && exact_back_jumps > 0
        ? array1d<dcomplex>(Nhs) : NULL,
             *sparse_next = sample_back_orbitals && exact_back_jumps > 0
        ? array1d<dcomplex>(Nhs) : NULL;
    vector<unsigned long long> sparse_marks(Nhs, 0);
    vector<int> sparse_active, sparse_next_active;
    sparse_active.reserve(Nhs);
    sparse_next_active.reserve(Nhs);
    unsigned long long sparse_generation = 0;
    if(sample_back_orbitals) {
      back_orbital_targets.resize(Nhs);
      back_orbital_factors.resize(Nhs);
      int *sample_state = array1d<int>(Nel);
      for(int s=0; s<Nhs; s++) {
        if(occ[s][0] == 0) {
          for(int k=0; k<Nvac; k++) {
            memcpy(sample_state, occ[s], Nel*sizeofint);
            const int orbital = virt[s][k];
            sample_state[0] = orbital;
            qsort(sample_state, Nel, sizeofint, intcmp);
            int position = 0;
            while(position < Nel && sample_state[position] != orbital) position++;
            int found = 0;
            const int target = binary_search3(sample_state, *occ, Nhs,
                Nel*sizeofint, &found, _mycmp3);
            if(!found || position >= Nel) {
              cerr<<"failed to build sampled backward 0->bath transition.\n";
              abort();
            }
            back_orbital_targets[s].push_back(target);
            back_orbital_factors[s].push_back(
                eo[position&1]*sqrtNvac*sqrt1_Nel);
          }
        } else {
          for(int k=0; k<Nel; k++) {
            memcpy(sample_state, occ[s], Nel*sizeofint);
            sample_state[k] = 0;
            qsort(sample_state, Nel, sizeofint, intcmp);
            int found = 0;
            const int target = binary_search3(sample_state, *occ, Nhs,
                Nel*sizeofint, &found, _mycmp3);
            if(!found) {
              cerr<<"failed to build sampled backward bath->0 transition.\n";
              abort();
            }
            back_orbital_targets[s].push_back(target);
            back_orbital_factors[s].push_back(
                eo[k&1]*sqrtNel*sqrt1_Nvac);
          }
        }
      }
      free1d(sample_state);
    }

    for(int n=0; n<ntraj_local; n++) {
      bzero(vec_for, vec_bytes);
      vec_for[initial_basis] = 1.0;
      dcomplex alp_for = alp_init,
               wgt_for = wgt_init;
      int excited_for = excited0,
          nj_for = 0;

      if(stratify_forward) {
        bzero(forward_jump_schedule, Jmax*sizeof(int));
        const double count_u = fmod(forward_count_shift +
            (trajectory_offset+n)*1.0/ntraj, 1.0);
        const int nj_forward = sepmb_sample_cdf(forward_count_cdf, count_u);
        if(sepmb_sample_jump_times_stratified(nstep, nj_forward,
              (unsigned long long)(trajectory_offset+n), forward_time_shift,
              stratify_single_jump, stratify_multi_jump,
              jumps_forward) != nj_forward) {
          cerr<<"failed to sample exact-orbital forward jump times.\n";
          abort();
        }
        for(int k=0; k<nj_forward; k++) {
          forward_jump_schedule[jumps_forward[k]] = 1;
        }
      }

      for(int j=1; j<=nstep; j++) {
        if(stratify_forward ? forward_jump_schedule[j]
                            : drand48() < jump_probability) {
          bzero(vec_tmp, vec_bytes);
          exc.multiply(vec_for, vec_tmp);
          for(int s=0; s<Nhs; s++) {
            vec_for[s] = vec_tmp[s]*inv_jump_normalization;
          }
          nj_for++;
          excited_for = 1-excited_for;
        }

        if(excited_for) {
          p0 = fpt*imag(alp_for);
          alp_for = (alp_for-dalp)*expfreqdt1 + dalp;
          pt = fpt*imag(alp_for);
          wgt_for *= exphwdt1*exp(0.5*I*(p0-pt)*delx);
        } else {
          wgt_for *= exphwdt1;
          alp_for *= expfreqdt1;
        }
        for(int s=0; s<Nhs; s++) vec_for[s] *= expEoccdt[1][s];
        jc[nj_for]++;

        const int iprb = j <= nwf ? measure_slot[j] : -1;
        if(iprb < 0) continue;

        const int parity = nj_for&1,
                  first_count = parity;
        const double back_accept_prb = back_accept[j][parity][first_count];
        if(back_accept_prb <= 0.0) continue;
        const double count_shift = drand48();

        for(int iback=0; iback<back_replicas; iback++) {
          const double count_u = fmod(count_shift +
              iback*inv_back_replicas, 1.0);
          const int nj_back = sepmb_sample_conditioned_count(
              back_accept[j][parity], j, parity, first_count, count_u);
          if(nj_back < 0) continue;
          const bool sparse_exact_back = sample_back_orbitals && exact_back_jumps > 0 &&
              nj_back <= exact_back_jumps;
          bzero(jumps_back, Jmax*sizeof(int));
          if(sepmb_sample_jump_times_stratified(j, nj_back,
                ((unsigned long long)(trajectory_offset+n)/nstep)*back_replicas+iback,
                backward_time_shift, stratify_single_jump,
                stratify_multi_jump,
                jumps_back) != nj_back) continue;

          sparse_active.clear();
          sparse_next_active.clear();
          if(sparse_exact_back) {
            sparse_back[initial_basis] = 1.0;
            sparse_active.push_back(initial_basis);
          } else if(!sample_back_orbitals) {
            bzero(vec_back, vec_bytes);
            vec_back[initial_basis] = 1.0;
          }
          dcomplex alp_back = alp_init,
                   wgt_back = wgt_init,
                   orbital_back_weight = 1.0;
          int excited_back = excited0,
              offset = 0,
              back_basis = initial_basis;

          for(int k=0; k<nj_back; k++) {
            const int nadvance = jumps_back[k]-1-offset;
            if(excited_back) {
              p0 = fpt*imag(alp_back);
              alp_back = (alp_back-dalp)*expfreqdt[nadvance] + dalp;
              pt = fpt*imag(alp_back);
              wgt_back *= exphwdt[nadvance]*exp(0.5*I*(p0-pt)*delx);
            } else {
              wgt_back *= exphwdt[nadvance];
              alp_back *= expfreqdt[nadvance];
            }
            if(sparse_exact_back) {
              for(int s : sparse_active) {
                sparse_back[s] *= expEoccdt[nadvance][s];
              }
              sparse_next_active.clear();
              sparse_generation++;
              for(int s : sparse_active) {
                const int degree = back_orbital_targets[s].size();
                for(int choice=0; choice<degree; choice++) {
                  const int target = back_orbital_targets[s][choice];
                  if(sparse_marks[target] != sparse_generation) {
                    sparse_marks[target] = sparse_generation;
                    sparse_next[target] = 0.0;
                    sparse_next_active.push_back(target);
                  }
                  sparse_next[target] += sparse_back[s]
                      *back_orbital_factors[s][choice]/(double)degree;
                }
              }
              for(int s : sparse_active) sparse_back[s] = 0.0;
              dcomplex *swap_buffer = sparse_back;
              sparse_back = sparse_next;
              sparse_next = swap_buffer;
              sparse_active.swap(sparse_next_active);
            } else if(sample_back_orbitals) {
              orbital_back_weight *= expEoccdt[nadvance][back_basis];
              const int degree = back_orbital_targets[back_basis].size();
              if(degree <= 0) {
                cerr<<"sampled backward determinant has no jump target.\n";
                abort();
              }
              int choice = (int)(drand48()*degree);
              if(choice >= degree) choice = degree-1;
              orbital_back_weight *= back_orbital_factors[back_basis][choice];
              back_basis = back_orbital_targets[back_basis][choice];
            } else {
              for(int s=0; s<Nhs; s++) {
                vec_back[s] *= expEoccdt[nadvance][s];
              }
              bzero(vec_tmp, vec_bytes);
              exc.multiply(vec_back, vec_tmp);
              for(int s=0; s<Nhs; s++) {
                vec_back[s] = vec_tmp[s]*inv_jump_normalization;
              }
            }
            offset = jumps_back[k]-1;
            excited_back = 1-excited_back;
          }

          const int nadvance = j-offset;
          if(excited_back) {
            p0 = fpt*imag(alp_back);
            alp_back = (alp_back-dalp)*expfreqdt[nadvance] + dalp;
            pt = fpt*imag(alp_back);
            wgt_back *= exphwdt[nadvance]*exp(0.5*I*(p0-pt)*delx);
          } else {
            wgt_back *= exphwdt[nadvance];
            alp_back *= expfreqdt[nadvance];
          }
          if(sparse_exact_back) {
            for(int s : sparse_active) {
              sparse_back[s] *= expEoccdt[nadvance][s];
            }
          } else if(sample_back_orbitals) {
            orbital_back_weight *= expEoccdt[nadvance][back_basis];
          } else {
            for(int s=0; s<Nhs; s++) {
              vec_back[s] *= expEoccdt[nadvance][s];
            }
          }

          if(excited_back != excited_for) {
            cerr<<"parity mismatch in exact-orbital SepMB projection.\n";
            abort();
          }
          const dcomplex ket_scale = wgt_for*Iton[nj_for%4]*sclf[j]
                                     *inv_ntraj*inv_back_replicas,
                         bra_scale = wgt_back*Iton[nj_back%4]
                                     *back_accept_prb*sclf[j];
          if(sparse_exact_back) {
            for(int s : sparse_active) {
              if(abs(vec_for[s]) <= 1.e-300 ||
                 abs(sparse_back[s]) <= 1.e-300) continue;
              csproj(ket_scale*vec_for[s], alp_for,
                     bra_scale*sparse_back[s], alp_back,
                     excited_for, basis_states[s], prb[iprb]);
            }
            for(int s : sparse_active) sparse_back[s] = 0.0;
          } else if(sample_back_orbitals) {
            if(abs(vec_for[back_basis]) > 1.e-300 &&
               abs(orbital_back_weight) > 1.e-300) {
              csproj(ket_scale*vec_for[back_basis], alp_for,
                     bra_scale*orbital_back_weight, alp_back,
                     excited_for, basis_states[back_basis], prb[iprb]);
            }
          } else {
            for(int s=0; s<Nhs; s++) {
              if(abs(vec_for[s]) <= 1.e-300 || abs(vec_back[s]) <= 1.e-300) continue;
              csproj(ket_scale*vec_for[s], alp_for,
                     bra_scale*vec_back[s], alp_back,
                     excited_for, basis_states[s], prb[iprb]);
            }
          }
        }
      }
    }

    free2d(expEoccdt);
    free1d(vec_for);
    if(vec_back) free1d(vec_back);
    if(sparse_back) free1d(sparse_back);
    if(sparse_next) free1d(sparse_next);
    free1d(vec_tmp);
    }
  } else {
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

    if(stratify_forward) {
      bzero(forward_jump_schedule, Jmax*sizeof(int));
      const double count_u = fmod(forward_count_shift +
          (trajectory_offset+n)*1.0/ntraj, 1.0);
      const int nj_forward = sepmb_sample_cdf(forward_count_cdf, count_u);
      if(sepmb_sample_jump_times_stratified(nstep, nj_forward,
            (unsigned long long)(trajectory_offset+n), forward_time_shift,
            stratify_single_jump, stratify_multi_jump,
            jumps_forward) != nj_forward) {
        cerr<<"failed to sample stratified forward jump times.\n";
        abort();
      }
      for(int k=0; k<nj_forward; k++) {
        forward_jump_schedule[jumps_forward[k]] = 1;
      }
    }

#ifdef _CHECK_PATH_
    path.clear();
#endif
    //the forward propagation
    for(int j=1; j<=nstep; j++) {
      //printf("    time: %4d", j);
      // jump seperatedly, switch the state between 0 and idx
      if(stratify_forward ? forward_jump_schedule[j]
                          : drand48() < jump_probability) {
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

      int iprb = j <= nwf ? measure_slot[j] : -1;
      if(iprb >= 0) {
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
        if(nj_min > j) continue;
        double back_accept_prb = back_accept[j][nj&1][nj_min];
        if(back_accept_prb <= 0.0) continue;
        dcomplex alp_for = alp,
                 wgt_for = wgt;
        set<int> state_for = state;
        const double count_shift = drand48();
        for(int iback=0; iback<back_replicas; iback++) {
#ifndef _CHECK_PATH_
          // Randomly shifted stratification keeps every replica marginally
          // uniform while spreading replicas over the conditional jump-count CDF.
          const double count_u = fmod(count_shift +
              iback*inv_back_replicas, 1.0);
          int nj_back = sepmb_sample_conditioned_count(
              back_accept[j][nj&1], j, nj&1, nj_min, count_u);
          if(nj_back < 0) continue;
          bzero(jumps_back, Jmax*sizeof(int));
          if(sepmb_sample_jump_times_stratified(j, nj_back,
                ((unsigned long long)(trajectory_offset+n)/nstep)*back_replicas+iback,
                backward_time_shift, stratify_single_jump,
                stratify_multi_jump,
                jumps_back) != nj_back) continue;

          int status = 0, fails = 0;
          do {
            status = (sampler.*PS[excited0][excited_for])(
                nj_back, S0, state_for, path);
            fails -= status;
          } while(status<0 && fails < 5);
          if(status<0) continue;
#else
          int nj_back = nj;
          cout<<"    "<<nj<<" jumps in the forward path"<<endl;
          print_set("     state : ", state_for, "\n");
          cout<<"    the backward path, t = "<<j<<endl;
#endif

          int offset = 0;
          alp     = alp_init;
          wgt     = wgt_init;
          excited = excited0;
          state   = S0;
          bool valid_path = ((int)path.size() == nj_back);
          for(int k=0; valid_path && k<nj_back; k++) {
            const int nadvance = jumps_back[k]-1-offset;
            if(excited) {
              p0   = fpt*imag(alp);
              alp  = (alp-dalp)*expfreqdt[nadvance] + dalp;
              pt   = fpt*imag(alp);
              wgt *= exphwdt[nadvance]*exp(0.5*I*(p0-pt)*delx);
            } else {
              wgt *= exphwdt[nadvance];
              alp *= expfreqdt[nadvance];
            }
            for(int l : state) wgt *= expEndt[nadvance][l];
            offset = jumps_back[k]-1;

            if(path[k].first < 0 || path[k].first >= Norb ||
               path[k].second < 0 || path[k].second >= Norb ||
               state.find(path[k].first) == state.end() ||
               state.find(path[k].second) != state.end()) {
              valid_path = false;
              break;
            }
            int pos = 0;
            if(excited) {
              state.erase(path[k].first);
              state.insert(path[k].second);
              idx = path[k].second;
              auto it = state.begin();
              for(; it != state.end() && *it != idx; ++it, ++pos) {}
              if(it == state.end() || (int)state.size() != Nel) {
                valid_path = false;
                break;
              }
            } else {
              idx = path[k].first;
              auto it = state.begin();
              for(; it != state.end() && *it != idx; ++it, ++pos) {}
              if(it == state.end()) {
                valid_path = false;
                break;
              }
              state.erase(path[k].first);
              state.insert(path[k].second);
              if((int)state.size() != Nel) {
                valid_path = false;
                break;
              }
            }
            sign = eo[pos%2];
#ifdef _TRACE_STATE_
            printf("  backward %2d-th jump for time %d, ", k, j);
            printf("     sign = %+d, switch %4d <-> %4d : \n",
                   sign, path[k].first, path[k].second);
#endif
            wgt *= sign*sqrtfct[excited];
            excited = 1-excited;
          }

          if(!valid_path) {
            alp = alp_for;
            wgt = wgt_for;
            state = state_for;
            excited = excited_for;
            continue;
          }

          if(excited) {
            p0   = fpt*imag(alp);
            alp  = (alp-dalp)*expfreqdt[j-offset] + dalp;
            pt   = fpt*imag(alp);
            wgt *= exphwdt[j-offset]*exp(0.5*I*(p0-pt)*delx);
          } else {
            wgt *= exphwdt[j-offset];
            alp *= expfreqdt[j-offset];
          }
          for(int l : state) wgt *= expEndt[j-offset][l];

          const double endpoint_prb = excited0
              ? sampler.get_Ptd(nj_back, d_kondo)
              : sampler.get_Qtd(nj_back, d_kondo);
          if(endpoint_prb <= 0.0) {
            alp = alp_for;
            wgt = wgt_for;
            state = state_for;
            excited = excited_for;
            continue;
          }
          const double measure = back_accept_prb*endpoint_prb*sclf[j];
          if(excited_for != excited) {
            cerr<<"parity mismatch between the forward and backward path.\n";
            abort();
          }
#ifdef _CHECK_PATH_
          print_set("     state : ", state, "\n");
          printf("wgt : %+1.16e %+1.16e %+1.16e %+1.16e\n",
                 real(wgt_for), imag(wgt_for), real(wgt), imag(wgt));
          printf("alp : %+1.16e %+1.16e %+1.16e %+1.16e\n",
                 real(alp_for), imag(alp_for), real(alp), imag(alp));
#endif
          csproj(wgt_for*Iton[nj%4]*sclf[j]*inv_ntraj*inv_back_replicas,
                 alp_for, wgt*Iton[nj_back%4]*measure, alp,
                 excited, state, prb[iprb]);

          alp = alp_for;
          wgt = wgt_for;
          state = state_for;
          excited = excited_for;
        }
      }
    }
	}
  }

#ifdef _YYY_MPI_
  double **avg  = array2d<double>(nmeas+1, Norb+3);
  MPI_Allreduce(*prb, *avg,  (1+nmeas)*(Norb+3),  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  memcpy(*prb, *avg, sizeof(double)*(1+nmeas)*(Norb+3));
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
    if(exact_orbitals && window_steps > 0) {
      fprintf(FL, "#PATCH_CHECK: SepMBpoisson v0.75 MPI-balanced phase-annihilated window propagation experimental\n");
    } else {
      fprintf(FL, exact_orbitals
          ? (sample_back_orbitals
              ? (exact_back_jumps > 0
                ? "#PATCH_CHECK: SepMBpoisson v0.71 direct-grid sampled-backward active\n"
                : "#PATCH_CHECK: SepMBpoisson v0.66 sampled-backward-orbital active\n")
              : (stratify_single_jump
                ? "#PATCH_CHECK: SepMBpoisson v0.64 exact-orbital 2D-single-jump-stratified active\n"
                : "#PATCH_CHECK: SepMBpoisson v0.62 exact-orbital Rao-Blackwell active\n"))
          : "#PATCH_CHECK: SepMBpoisson v0.61 exact-bernoulli stratified-backward active\n");
    }
    fprintf(FL, "#discretizing the bath:\n");
    for(int n=0; n<Norb; n++) fprintf(FL, "#%6d %1.16e %1.16e\n", n, cpl[n], En[n]);
    fprintf(FL, "#sampling: jump_strength=%1.16e jump_probability=%1.16e log_scale=%1.16e back_replicas=%d stratify_forward=%d exact_orbitals=%d stratify_single_jump_time=%d stratify_multi_jump_time=%d sample_back_orbitals=%d exact_back_jumps=%d\n",
            jump_strength, jump_probability, log_scale, back_replicas,
            stratify_forward ? 1 : 0, exact_orbitals ? 1 : 0,
            stratify_single_jump ? 1 : 0, stratify_multi_jump ? 1 : 0,
            sample_back_orbitals ? 1 : 0,
            exact_back_jumps);
    fprintf(FL, "#adaptive measurement: gap=%1.16e period_steps=%d nmeas=%d nwf=%d tmax=%1.16e forced_stride=%d window_steps=%d window_count=%d window_pair_replicas=%d window_alpha_bin=%1.16e\n",
            gap, period_steps, nmeas, nwf, nwf*dt, forced_measure_stride,
            window_steps, window_steps > 0 ? (nwf+window_steps-1)/window_steps : 0,
            window_pair_replicas, window_alpha_bin);
    double *rlt = array1d<double>(Norb+3);
    int il = 0;
    for(int t=0; t<=nwf; t++)  {
      while(il+1 <= nmeas && measure_steps[il+1] < t) il++;
      int ir = il+1 <= nmeas ? il+1 : il,
          tl = measure_steps[il],
          tr = measure_steps[ir];
      if(il == ir || tr <= tl) {
        for(int k=0; k<Norb+3; k++) rlt[k] = prb[il][k];
      } else {
        double f = (t-tl)*1.0/(tr-tl);
        for(int k=0; k<Norb+3; k++) rlt[k] = (1.0-f)*prb[il][k] + f*prb[ir][k];
      }

      double norm = rlt[0]/Nel;
      if(fabs(norm) > 1.e-300) {
        fprintf(FL, "%12.8f %+1.16e %+1.16e %+1.16e", t*dt, (double)Nel, rlt[1]/norm, rlt[2]/norm);
        for(int k=3; k<Norb+3; k++) fprintf(FL, " %+1.16e", rlt[k]/norm);
      } else {
        fprintf(FL, "%12.8f %+1.16e %+1.16e %+1.16e", t*dt, rlt[0], rlt[1], rlt[2]);
        for(int k=3; k<Norb+3; k++) fprintf(FL, " %+1.16e", rlt[k]);
      }
      fprintf(FL, "\n");
    }
    free1d(rlt);
    fclose(FL);

    sprintf(fnm, "ahm-jcmb-s%d-n%d-%d.dat", Norb, Nel, ntraj);
    FL = fopen(fnm, "w");
    for(int t=0; t<nstep; t++) {
      fprintf(FL, "%20d %1.16e\n", t, jc[t]*1.0/ntraj);
    }
    fclose(FL);
  }

  free3d(back_accept);
  free2d(prb);
  free1d(measure_slot);
  free1d(sclf);
  free1d(jumps_back);
  free1d(jumps_forward);
  free1d(forward_jump_schedule);
  free1d(jc);

	return;
}
