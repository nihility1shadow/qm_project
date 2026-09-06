#!/bin/bash
set -euo pipefail
project=/data/home/yd101802/yd101802/nonadia
cd "$project"
manifest="$project/cloud-runs/v113-4000-20260906.tsv"
if [[ -e "$manifest" ]]; then
  cat "$manifest"
  exit 1
fi
printf 'role\tjob\treference_job\tseed\tpaths\tfock_levels\tdistance\tranks\n' > "$manifest"
ref384=$(sbatch --parsable -J v113_ref4000_f384 --nodes=2 --ntasks=128 --ntasks-per-node=64 --mem=24G --time=48:00:00 --export=ALL,BINARY=na_mpi_v113_cache.out,FLEVELS=384 run_poisson_lowmem.slurm 8000 100 30 15 3 800000 1 1134000)
printf 'reference384\t%s\t-\t1134000\t100\t384\t3\t128\n' "$ref384" >> "$manifest"
ref512=$(sbatch --parsable -J v113_ref4000_f512 --nodes=2 --ntasks=128 --ntasks-per-node=64 --mem=24G --time=48:00:00 --export=ALL,BINARY=na_mpi_v113_cache.out,FLEVELS=512 run_poisson_lowmem.slurm 8000 100 30 15 3 800000 1 1134000)
printf 'reference512\t%s\t-\t1134000\t100\t512\t3\t128\n' "$ref512" >> "$manifest"
for seed in 1134001 1134002 1134003; do
  cache="$project/cloud-runs/v111-distributed-20260905/$ref384/reference-observables.dat"
  child=$(sbatch --parsable -J "v113_4000_s$seed" --nodes=2 --ntasks=128 --ntasks-per-node=64 --mem=4G --time=48:00:00 --dependency="afterok:$ref384" --export="ALL,BINARY=na_mpi_v113_cache.out,FLEVELS=384,SEP_MB_REFERENCE_CACHE=$cache" run_poisson_lowmem.slurm 8000 1000000 30 15 3 800000 1 "$seed")
  printf 'sample384\t%s\t%s\t%s\t1000000\t384\t3\t128\n' "$child" "$ref384" "$seed" >> "$manifest"
done
cache="$project/cloud-runs/v111-distributed-20260905/$ref512/reference-observables.dat"
child=$(sbatch --parsable -J v113_4000_f512 --nodes=2 --ntasks=128 --ntasks-per-node=64 --mem=4G --time=48:00:00 --dependency="afterok:$ref512" --export="ALL,BINARY=na_mpi_v113_cache.out,FLEVELS=512,SEP_MB_REFERENCE_CACHE=$cache" run_poisson_lowmem.slurm 8000 1000000 30 15 3 800000 1 1134001)
printf 'sample512\t%s\t%s\t1134001\t1000000\t512\t3\t128\n' "$child" "$ref512" >> "$manifest"
refd2=$(sbatch --parsable -J v113_ref4000_d2 --nodes=2 --ntasks=128 --ntasks-per-node=64 --mem=4G --time=12:00:00 --export=ALL,BINARY=na_mpi_v113_cache.out,FLEVELS=384 run_poisson_lowmem.slurm 8000 100 30 15 2 65536 1 1134000)
printf 'reference_d2\t%s\t-\t1134000\t100\t384\t2\t128\n' "$refd2" >> "$manifest"
cat "$manifest"
