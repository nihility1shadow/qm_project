#!/bin/bash
set -euo pipefail
cd /data/home/yd101802/yd101802/nonadia
manifest=cloud-runs/v114-signed-csr-short-20260906.tsv
test ! -e "$manifest" || { echo 'Existing manifest; inspect jobs instead of duplicating.'; exit 2; }
printf 'role\tjob\n' > "$manifest"
submit() {
  local role="$1"; shift
  local id
  id=$(sbatch --parsable "$@")
  printf '%s\t%s\n' "$role" "$id" | tee -a "$manifest"
}
submit small --partition=cmh --nodes=1 --ntasks=4 --time=00:20:00 --mem=1G --export=ALL,BINARY=na_mpi_v114_signed_csr.out,FLEVELS=384,SEP_MB_REFERENCE_SIGNED_CSR=1 run_poisson_lowmem.slurm 200 1000 10 5 3 800000 1 111001
submit small_fallback --partition=cmh --nodes=1 --ntasks=4 --time=00:20:00 --mem=1G --export=ALL,BINARY=na_mpi_v114_signed_csr.out,FLEVELS=384,SEP_MB_REFERENCE_SIGNED_CSR=0 run_poisson_lowmem.slurm 200 1000 10 5 3 800000 1 111001
submit D2_64 --partition=cmh --nodes=1 --ntasks=64 --time=00:20:00 --mem=4G --export=ALL,BINARY=na_mpi_v114_signed_csr.out,FLEVELS=384,SEP_MB_REFERENCE_SIGNED_CSR=1 run_poisson_lowmem.slurm 20 100 30 15 2 800000 1 111001
submit D3_64 --partition=cmh --nodes=1 --ntasks=64 --time=00:20:00 --mem=24G --export=ALL,BINARY=na_mpi_v114_signed_csr.out,FLEVELS=384,SEP_MB_REFERENCE_SIGNED_CSR=1 run_poisson_lowmem.slurm 20 100 30 15 3 800000 1 111001
submit D3_128 --partition=cmh --nodes=2 --ntasks=128 --time=00:20:00 --mem=24G --export=ALL,BINARY=na_mpi_v114_signed_csr.out,FLEVELS=384,SEP_MB_REFERENCE_SIGNED_CSR=1 run_poisson_lowmem.slurm 20 100 30 15 3 800000 1 111001
submit cache_read --partition=cmh --nodes=1 --ntasks=4 --time=00:20:00 --mem=1G --export=ALL,BINARY=na_mpi_v114_signed_csr.out,FLEVELS=384,SEP_MB_REFERENCE_SIGNED_CSR=1,SEP_MB_REFERENCE_CACHE=/data/home/yd101802/yd101802/nonadia/cloud-runs/v111-distributed-20260905/647662/reference-observables.dat run_poisson_lowmem.slurm 200 1000 10 5 3 800000 1 111001
