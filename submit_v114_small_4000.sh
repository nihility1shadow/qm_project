#!/bin/bash
set -euo pipefail
cd /data/home/yd101802/yd101802/nonadia
manifest=cloud-runs/v114-small-4000-20260906.tsv
test ! -e "$manifest" || { echo 'Manifest already exists; inspect current jobs.'; exit 2; }
printf 'role\tjob\treference_job\tseed\tpaths\tranks\n' > "$manifest"
common=(--partition=cmh --nodes=1 --ntasks=64 --ntasks-per-node=64 --time=12:00:00 --mem=2G)
producer=$(sbatch --parsable "${common[@]}" --export=ALL,BINARY=na_mpi_v114_signed_csr.out,FLEVELS=384,SEP_MB_REFERENCE_SIGNED_CSR=1 run_poisson_lowmem.slurm 8000 100 10 5 3 4096 1 1144000)
printf 'reference\t%s\t-\t1144000\t100\t64\n' "$producer" | tee -a "$manifest"
cache="/data/home/yd101802/yd101802/nonadia/cloud-runs/v111-distributed-20260905/$producer/reference-observables.dat"
for seed in 1144001 1144002 1144003; do
  job=$(sbatch --parsable "${common[@]}" --dependency="afterok:$producer" --export="ALL,BINARY=na_mpi_v114_signed_csr.out,FLEVELS=384,SEP_MB_REFERENCE_SIGNED_CSR=1,SEP_MB_REFERENCE_CACHE=$cache" run_poisson_lowmem.slurm 8000 1000000 10 5 3 4096 1 "$seed")
  printf 'sample\t%s\t%s\t%s\t1000000\t64\n' "$job" "$producer" "$seed" | tee -a "$manifest"
done
