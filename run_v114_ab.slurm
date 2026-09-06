#!/bin/bash
#SBATCH --partition=cmh
#SBATCH --nodes=2
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=64
#SBATCH --time=00:30:00
#SBATCH --mem=24G
set -euo pipefail
cd /data/home/yd101802/yd101802/nonadia
export FLEVELS=384 SEP_MB_REFERENCE_SIGNED_CSR=1
for variant in old_first new_first new_second old_second; do
  export BENCH_CASE="$variant"
  case "$variant" in
    old*) export BINARY=na_mpi_v113_cache.out;;
    new*) export BINARY=na_mpi_v114_signed_csr.out;;
  esac
  bash run_v114_ab_case.sh 20 100 30 15 3 800000 1 111001
done
