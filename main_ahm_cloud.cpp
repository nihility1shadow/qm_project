#include <mpi.h>
#include <sys/time.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <algorithm>
#include <climits>
#include <cstring>
#include <sys/resource.h>

#include "ahm.h"
#include "myid.h"

int myid = 0;
int master = 0;
int nproc = 1;

namespace {
const double EV_PER_HARTREE = 27.211386245988;

unsigned long mix_seed(unsigned long value) {
  value ^= value >> 16;
  value *= 0x7feb352dUL;
  value ^= value >> 15;
  value *= 0x846ca68bUL;
  value ^= value >> 16;
  return value;
}

void fail(const char *message) {
  if (myid == master) std::fprintf(stderr, "%s\n", message);
  MPI_Abort(MPI_COMM_WORLD, 1);
}
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  master = nproc - 1;

  if (argc != 6 || std::strtol(argv[1], NULL, 10) != 10) {
    fail("Syntax: na_mpi_cloud.out 10 job ntraj Norb Nel");
  }

  const int job = std::strtol(argv[2], NULL, 10);
  const int ntraj = std::strtol(argv[3], NULL, 10);
  const int Norb = std::strtol(argv[4], NULL, 10);
  const int Nel = std::strtol(argv[5], NULL, 10);
  if (job != 0 && job != 1) fail("The cloud AHM target supports job 0 or 1.");
  if (job == 0 && nproc != 1) {
    fail("Job 0 writes a single QM trajectory and must use one MPI task.");
  }
  if (ntraj <= 0 || Norb < 3 || Nel <= 0 || Nel >= Norb) {
    fail("Need ntraj > 0, Norb >= 3, and 0 < Nel < Norb.");
  }

  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long seed_base = 0;
  if (myid == master) {
    const char *env_seed = std::getenv("AHM_SEED");
    if (env_seed) {
      seed_base = std::strtoul(env_seed, NULL, 10);
    } else {
      seed_base = static_cast<unsigned long>(tv.tv_sec);
      seed_base ^= static_cast<unsigned long>(tv.tv_usec) << 11;
      const char *slurm_job_id = std::getenv("SLURM_JOB_ID");
      if (slurm_job_id) {
        seed_base ^= std::strtoul(slurm_job_id, NULL, 10) * 0x9e3779b1UL;
      }
      seed_base = mix_seed(seed_base);
    }
  }
  MPI_Bcast(&seed_base, 1, MPI_UNSIGNED_LONG, master, MPI_COMM_WORLD);
  const unsigned long rank_seed =
      mix_seed(seed_base ^ (static_cast<unsigned long>(myid) + 1) * 0x85ebca6bUL);
  srand48(static_cast<long>(rank_seed & 0x7fffffffUL));

  double wc = 10.0 / EV_PER_HARTREE;
  double eta = 0.01;
  double delE = -0.6319;
  const double delx = 2.0;
  const double xmin = -8.0;
  const double xmax = 8.0;
  const double pinit = 0.0;
  const double mass = 14583.1067146087;
  const double freq = 3.6749323758566211e-03;
  const int npt = 1024;
  int nstep = 8000;
  double dt = 0.5;

  const char *env_wc_ev = std::getenv("AHM_WC_EV");
  const char *env_eta = std::getenv("AHM_ETA");
  const char *env_dele_ev = std::getenv("AHM_DELE_EV");
  const char *env_nstep = std::getenv("AHM_NSTEP");
  const char *env_bath_model = std::getenv("AHM_BATH_MODEL");
  const char *bath_model = env_bath_model ? env_bath_model : "oxygen";
  if (env_wc_ev) wc = std::atof(env_wc_ev) / EV_PER_HARTREE;
  if (env_eta) eta = std::atof(env_eta);
  if (env_dele_ev) delE = std::atof(env_dele_ev) / EV_PER_HARTREE;
  if (env_nstep) nstep = std::atoi(env_nstep);
  const char *env_dt = std::getenv("AHM_DT");
  if (env_dt) dt = std::atof(env_dt);
  if (wc <= 0.0 || wc > 10.0 / EV_PER_HARTREE) {
    fail("AHM_WC_EV must be in (0, 10] eV for this scan.");
  }
  if (nstep <= 0) fail("AHM_NSTEP must be positive.");
  if (!std::isfinite(dt) || dt <= 0.0) fail("AHM_DT must be finite and positive.");
  if (!std::isfinite(eta) || eta < 0.0 || !std::isfinite(delE) || !std::isfinite(wc))
    fail("Physical energies and coupling parameters must be finite; eta must be nonnegative.");
  if (std::strcmp(bath_model, "oxygen") != 0 &&
      std::strcmp(bath_model, "semicircle") != 0)
    fail("AHM_BATH_MODEL must be either oxygen or semicircle.");
  const char *env_path_local_basis = std::getenv("SEP_MB_PATH_LOCAL_BASIS");
  const bool path_local_basis = job == 1 && env_path_local_basis &&
      std::atoi(env_path_local_basis) != 0;
  long double hilbert = 1.0L, legacy_numerator = 1.0L;
  for (int k=1; k<=std::min(Nel,Norb-Nel); ++k)
    hilbert *= (long double)(Norb-k+1)/k;
  for (int k=1; k<=Nel; ++k) legacy_numerator *= Norb-k+1;
  if (myid == master) {
    std::printf("#AHAU_PREFLIGHT orbitals=%d electrons=%d determinants_estimate=%.0Lf "
                "grid_points=%d single_full_wavefunction_gib=%.6Lf "
                "full_basis_requested=%d\n", Norb, Nel, hilbert, npt,
                hilbert*npt*sizeof(dcomplex)/(1024.0L*1024*1024), path_local_basis?0:1);
    std::fflush(stdout);
  }
  const char *env_preflight = std::getenv("AHM_PREFLIGHT_ONLY");
  if (env_preflight && std::atoi(env_preflight)) { MPI_Finalize(); return 0; }
  if (!path_local_basis &&
      (hilbert*npt > INT_MAX || legacy_numerator > INT_MAX)) {
    fail("Full-basis request exceeds the legacy integer-safe layout. "
         "Use job 1 with SEP_MB_PATH_LOCAL_BASIS=1; this is not a completed full-space QM run.");
  }
  const double computation_started = MPI_Wtime();

  if (myid == master) {
    std::printf("#AHAU_MPI nproc=%d seed_base=%lu bath_model=%s wc_eV=%1.16e "
                "eta_au2=%1.16e eta_eV2=%1.16e delE_eV=%1.16e nstep=%d dt=%g\n",
                nproc, seed_base,
                bath_model, wc * EV_PER_HARTREE, eta,
                eta * EV_PER_HARTREE * EV_PER_HARTREE,
                delE * EV_PER_HARTREE, nstep, dt);
  }

  AHM ahm;
  ahm.set_mass(mass);
  ahm.set_freq(freq);
  ahm.set_delx(delx);
  ahm.set_grids(npt, xmin, xmax);
  if (std::strcmp(bath_model, "semicircle") == 0)
    ahm.dissemicircle(Norb, eta, wc);
  else
    ahm.diseven(Norb, eta, wc);
  ahm.set_Nel(Nel, Norb);
  ahm.set_delE(delE);
  if (!path_local_basis) {
    ahm.set_basis();
    ahm.calc_Eocc();
    ahm.set_exc();
    ahm.set_Nex();
  } else if (myid == master) {
    std::printf("#AHAU_BASIS mode=path-local full_determinant_basis=skipped\n");
  }

  const double xinit = delx;
  const dcomplex alp0 = std::sqrt(mass * freq / 2.0) *
                        (xinit + (pinit / (mass * freq)) * I);
  if (job == 0) {
    ahm.qm(nstep, dt, alp0);
  } else {
    std::set<int> state;
    for (int j = 0; j < Nel; ++j) state.insert(j);
    ahm.SepMBpoisson(ntraj, nstep, dt, alp0, 1.0, state);
  }

  struct rusage usage;
  const bool usage_ok = getrusage(RUSAGE_SELF, &usage) == 0;
  long peak = usage_ok ? usage.ru_maxrss : -1, maximum_peak=0, peak_sum=0;
  double cpu[2] = {usage_ok ? usage.ru_utime.tv_sec+usage.ru_utime.tv_usec*1.e-6 : -1.0,
                   usage_ok ? usage.ru_stime.tv_sec+usage.ru_stime.tv_usec*1.e-6 : -1.0}, total_cpu[2];
  const double elapsed = MPI_Wtime()-computation_started;
  double maximum_elapsed=0.0;
  MPI_Reduce(&peak,&maximum_peak,1,MPI_LONG,MPI_MAX,master,MPI_COMM_WORLD);
  MPI_Reduce(&peak,&peak_sum,1,MPI_LONG,MPI_SUM,master,MPI_COMM_WORLD);
  MPI_Reduce(cpu,total_cpu,2,MPI_DOUBLE,MPI_SUM,master,MPI_COMM_WORLD);
  MPI_Reduce(&elapsed,&maximum_elapsed,1,MPI_DOUBLE,MPI_MAX,master,MPI_COMM_WORLD);
  if(myid==master) std::printf("#AHAU_RESOURCES ranks=%d max_rank_peak_rss_kb=%ld "
      "sum_rank_peak_rss_kb=%ld cpu_user_seconds_sum=%.6f cpu_system_seconds_sum=%.6f "
      "computation_wall_seconds_max=%.6f (peak sum is not simultaneous RSS)\n",
      nproc,maximum_peak,peak_sum,total_cpu[0],total_cpu[1],maximum_elapsed);
  MPI_Finalize();
  return 0;
}
