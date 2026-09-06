#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=4G
set -euo pipefail
module purge
module load gcc/12.1.0 openmpi/4.1.8
project=/data/home/yd101802/yd101802/nonadia
run="$project/cloud-runs/v111-distributed-20260905/$SLURM_JOB_ID"
mkdir -p "$run"
cd "$run"
export LD_LIBRARY_PATH="$project/.deps/install/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export AHM_WC_EV=2.7 AHM_ETA=0.30e-6 AHM_DELE_EV=-17.194874968839816
export AHM_NSTEP=8000 AHM_SEED=1134001
binary="$project/na_mpi_v113_csr.out"
env | sort | grep -E '^(AHM_|SLURM_NTASKS)' > config.txt
sha256sum "$binary" > binary.sha256
/usr/bin/time -f 'wall_seconds=%e launcher_peak_rss_kb=%M' -o time.txt \
  mpirun --bind-to core "$binary" 10 0 1 10 5 > program.out 2> program.err
