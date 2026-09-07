set -euo pipefail
cd /data/home/yd101802/yd101802/nonadia
manifest=cloud-runs/sampling-comparison-4000-20260907.tsv
test ! -e "$manifest" || { cat "$manifest"; exit 2; }
states=$(squeue -r -h -j 647806 -o %T | sort -u)
if [ "$states" != PENDING ]; then
  echo "Old array is not entirely pending; inspect before replacing it."
  exit 3
fi
scancel 647806
wrap=$(cat <<'RUN_SCRIPT'
set -eu
module purge
module load gcc/12.1.0 openmpi/4.1.8
project=/data/home/yd101802/yd101802/nonadia
index=$SLURM_ARRAY_TASK_ID
case "$index" in 0|1) rate=1.5; stepmode=1;; 2|3) rate=1.0; stepmode=0;; 4|5) rate=1.5; stepmode=0;; *) exit 2;; esac
run=$project/cloud-runs/sampling-comparison-4000-20260907/${SLURM_ARRAY_JOB_ID}_${index}
mkdir -p "$run"
cd "$run"
export LD_LIBRARY_PATH="$project/.deps/install/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export AHM_WC_EV=2.7 AHM_ETA=0.30e-6 AHM_DELE_EV=-17.194874968839816
export AHM_NSTEP=8000 AHM_SEED=$((1159001 + index % 2))
export SEP_MB_TMAX=4000 SEP_MB_DETERMINISTIC_FOCK=0 SEP_MB_EXACT_ORBITALS=0
export SEP_MB_PATH_LOCAL_BASIS=1 SEP_MB_REFERENCE_OUTPUT=1 SEP_MB_REFERENCE_GRID_SPLIT=1
export SEP_MB_REFERENCE_SPLIT_OPERATOR=1 SEP_MB_REFERENCE_DISTANCE=3 SEP_MB_REFERENCE_MAX_STATES=800000
export SEP_MB_REFERENCE_FOCK_STATES=384 SEP_MB_REFERENCE_MPI=1 SEP_MB_REFERENCE_SIGNED_CSR=1
export SEP_MB_REFERENCE_CACHE=$project/cloud-runs/v111-distributed-20260905/647702/reference-observables.dat
export SEP_MB_RATE_SCALE="$rate" SEP_MB_MEASURE_STRIDE=4
export SEP_MB_BACK_REPLICAS=16 SEP_MB_BACK_REPLICAS_MIN=16 SEP_MB_BACK_REPLICA_POWER=2
export SEP_MB_STRATIFY_FORWARD_COUNT=1 SEP_MB_STRATIFY_SINGLE_JUMP_TIME=1 SEP_MB_SAMPLE_BACK_ORBITALS=1
export SEP_MB_STRATIFY_FORWARD_ORBITALS=1 SEP_MB_STRATIFY_FORWARD_STEPS="$stepmode" SEP_MB_STRATIFY_BACK_PATHS=1
export SEP_MB_RQMC_REPLICATES=4 SEP_MB_EXACT_BACK_JUMPS=0 SEP_MB_ALL_ORDER_BACK_DP=1
export SEP_MB_RECURRENCE_DENSE=1 SEP_MB_RECURRENCE_STRIDE=3
binary=$project/na_mpi_v114_signed_csr.out
env | sort | grep -E '^(AHM_|SEP_MB_|SLURM_NTASKS)' > config.txt
sha256sum "$binary" > binary.sha256
/usr/bin/time -f 'wall_seconds=%e launcher_peak_rss_kb=%M' -o time.txt mpirun --bind-to core --map-by ppr:64:node bash -c 'exec /usr/bin/time -f "%U %S %e %M" -o "process-${OMPI_COMM_WORLD_RANK}.txt" "$@"' _ "$binary" 10 1 100000 30 15 > program.out 2> program.err
awk '{u+=$1;s+=$2;if($3>w)w=$3;if($4>m)m=$4;r+=$4;n++} END{printf "ranks=%d user_cpu_seconds_sum=%.6f system_cpu_seconds_sum=%.6f max_process_wall_seconds=%.6f max_rank_peak_rss_kb=%.0f sum_rank_peak_rss_kb=%.0f\n",n,u,s,w,m,r}' process-[0-9]*.txt > process-resources.txt
RUN_SCRIPT
)
job=$(sbatch --parsable --partition=cmh --nodes=2 --ntasks=128 --ntasks-per-node=64 --array=0-5%1 --time=02:00:00 --mem=4G --job-name=sampling_30e15_4000 --wrap="$wrap")
printf 'array_job\ttask\tforward_steps\trate_scale\tseed\tpaths\tranks\treference_job\n' > "$manifest"
for index in 0 1 2 3 4 5; do
  case "$index" in 0|1) rate=1.5; stepmode=1;; 2|3) rate=1.0; stepmode=0;; 4|5) rate=1.5; stepmode=0;; esac
  printf '%s\t%s\t%s\t%s\t%s\t100000\t128\t647702\n' "$job" "$index" "$stepmode" "$rate" "$((1159001+index%2))" | tee -a "$manifest"
done