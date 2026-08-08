#include <mpi.h>
#include <sys/time.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>

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
  const double dt = 0.5;

  const char *env_wc_ev = std::getenv("AHM_WC_EV");
  const char *env_eta = std::getenv("AHM_ETA");
  const char *env_dele_ev = std::getenv("AHM_DELE_EV");
  const char *env_nstep = std::getenv("AHM_NSTEP");
  if (env_wc_ev) wc = std::atof(env_wc_ev) / EV_PER_HARTREE;
  if (env_eta) eta = std::atof(env_eta);
  // Preserve the current server model's 0.6 conversion factor exactly.
  if (env_dele_ev) delE = std::atof(env_dele_ev) * 0.6 / EV_PER_HARTREE;
  if (env_nstep) nstep = std::atoi(env_nstep);
  if (wc <= 0.0 || wc > 10.0 / EV_PER_HARTREE) {
    fail("AHM_WC_EV must be in (0, 10] eV for this scan.");
  }
  if (nstep <= 0) fail("AHM_NSTEP must be positive.");

  if (myid == master) {
    std::printf("#AHAU_MPI nproc=%d seed_base=%lu wc_eV=%1.16e eta=%1.16e "
                "delE_eV=%1.16e nstep=%d dt=%g\n",
                nproc, seed_base,
                wc * EV_PER_HARTREE, eta,
                delE * EV_PER_HARTREE, nstep, dt);
  }

  AHM ahm;
  ahm.set_mass(mass);
  ahm.set_freq(freq);
  ahm.set_delx(delx);
  ahm.set_grids(npt, xmin, xmax);
  ahm.diseven(Norb, eta, wc);
  ahm.set_Nel(Nel, Norb);
  ahm.set_delE(delE);
  ahm.set_basis();
  ahm.calc_Eocc();
  ahm.set_exc();
  ahm.set_Nex();

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

  MPI_Finalize();
  return 0;
}
