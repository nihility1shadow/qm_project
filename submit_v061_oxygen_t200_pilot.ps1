$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = "E:\project_qm"
$ScanRoot = Join-Path $ProjectRoot "scan_v061_oxygen_t200_q10_20260810"
$CandidateFile = Join-Path $ScanRoot "pilot_candidates.csv"
$SubmissionFile = Join-Path $ScanRoot "pilot_submission.csv"
$RemoteRoot = "/data/home/yd101802/yd101802/nonadia"
$RemoteScanRoot = "$RemoteRoot/cloud-runs/v061-oxygen-t200-q10-20260810/pilot"
$Key = Join-Path $env:TEMP "qm_codex_ssh\qm_codex_ed25519"
$KnownHosts = Join-Path $ProjectRoot "tmp\cloud_known_hosts_20260809"
$HostName = "yd101802@36.212.18.57"

$DelE = -17.194874968839816
$Ntraj = 8000000
$Nstep = 400
$Ntasks = 64
$Repeats = 3
$RandomSeed = 20260810

New-Item -ItemType Directory -Force -Path $ScanRoot | Out-Null

# The center point measures the variance reduction from extra backward replicas.
# Random probes are reproducible and expand the parameter search each round.
$candidates = @(
    [pscustomobject]@{ case_id = "p01_center_b04"; param_id = "center"; origin = "replica_control"; random_seed = 0; wc_eV = 4.0;      eta = 0.000170000; back_replicas = 4  },
    [pscustomobject]@{ case_id = "p02_center_b08"; param_id = "center"; origin = "replica_control"; random_seed = 0; wc_eV = 4.0;      eta = 0.000170000; back_replicas = 8  },
    [pscustomobject]@{ case_id = "p03_center_b12"; param_id = "center"; origin = "replica_control"; random_seed = 0; wc_eV = 4.0;      eta = 0.000170000; back_replicas = 12 },
    [pscustomobject]@{ case_id = "p04_center_b16"; param_id = "center"; origin = "replica_control"; random_seed = 0; wc_eV = 4.0;      eta = 0.000170000; back_replicas = 16 },
    [pscustomobject]@{ case_id = "p05_random_01"; param_id = "random01"; origin = "random"; random_seed = $RandomSeed; wc_eV = 4.145914; eta = 0.000124380; back_replicas = 12 },
    [pscustomobject]@{ case_id = "p06_random_02"; param_id = "random02"; origin = "random"; random_seed = $RandomSeed; wc_eV = 3.509066; eta = 0.000124874; back_replicas = 12 },
    [pscustomobject]@{ case_id = "p07_random_03"; param_id = "random03"; origin = "random"; random_seed = $RandomSeed; wc_eV = 4.689026; eta = 0.000207282; back_replicas = 12 }
)
$candidates | Export-Csv -NoTypeInformation -Encoding utf8 $CandidateFile

$rows = [System.Collections.Generic.List[object]]::new()
if (Test-Path $SubmissionFile) {
    foreach ($row in Import-Csv $SubmissionFile) { $rows.Add($row) }
}

function Invoke-Remote([string]$Command) {
    for ($attempt = 1; $attempt -le 30; $attempt++) {
        $previousPreference = $ErrorActionPreference
        $ErrorActionPreference = "Continue"
        $result = & ssh -i $Key -o BatchMode=yes -o UserKnownHostsFile=$KnownHosts `
            -p 22 $HostName $Command 2>&1
        $sshExitCode = $LASTEXITCODE
        $ErrorActionPreference = $previousPreference
        if ($sshExitCode -eq 0) {
            $jobId = $result | ForEach-Object { "$_" } |
                Where-Object { $_ -match '^\d+$' } | Select-Object -Last 1
            if ($jobId) { return $jobId.Trim() }
        }
        if ($attempt -lt 30) { Start-Sleep -Seconds 5 }
    }
    throw "Submission failed after 30 attempts"
}

function Save-Row([object]$Row) {
    $script:rows.Add($Row)
    $script:rows | Export-Csv -NoTypeInformation -Encoding utf8 $SubmissionFile
}

function Already-Submitted([string]$CaseId, [string]$Method, [int]$Repeat) {
    return [bool]($rows | Where-Object {
        $_.case_id -eq $CaseId -and $_.method -eq $Method -and
        [int]$_.repeat -eq $Repeat
    })
}

$candidateNumber = 0
foreach ($candidate in $candidates) {
    $candidateNumber++
    $caseId = $candidate.case_id
    $wc = [double]$candidate.wc_eV
    $eta = [double]$candidate.eta
    $backReplicas = [int]$candidate.back_replicas

    foreach ($repeat in 1..$Repeats) {
        if (Already-Submitted $caseId "v061" $repeat) { continue }
        $runDir = "$RemoteScanRoot/$caseId/rep$repeat"
        $jobName = "v61p{0:d2}r$repeat" -f $candidateNumber
        $wrap = "set -euo pipefail; module purge; module load gcc/12.1.0; " +
            "module load openmpi/4.1.8; export OMP_NUM_THREADS=1 " +
            "MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; " +
            "export LD_LIBRARY_PATH=$RemoteRoot/.deps/install/lib:" +
            '$LD_LIBRARY_PATH' + "; mpirun --bind-to core --map-by core " +
            "-np $Ntasks $RemoteRoot/na_mpi_v061.out 10 1 $Ntraj 10 5 " +
            "> program.out 2>&1"
        $remote = "cd $RemoteRoot && mkdir -p $runDir && " +
            "sbatch --parsable --partition=cmh --ntasks=$Ntasks " +
            "--cpus-per-task=1 --mem-per-cpu=1G --time=08:00:00 " +
            "--job-name=$jobName --chdir=$runDir " +
            "--export=ALL,AHM_WC_EV=$wc,AHM_ETA=$eta,AHM_DELE_EV=$DelE," +
            "AHM_NSTEP=$Nstep,SEP_MB_BACK_REPLICAS=$backReplicas," +
            "SEP_MB_STRATIFY_FORWARD_COUNT=0 " +
            "--output=slurm-%j.out --error=slurm-%j.err --wrap='$wrap'"
        $jobId = Invoke-Remote $remote
        Save-Row ([pscustomobject]@{
            stage = "pilot"; case_id = $caseId; param_id = $candidate.param_id
            origin = $candidate.origin; random_seed = $candidate.random_seed
            method = "v061"; repeat = $repeat; wc_eV = $wc; eta = $eta
            delE_eV = $DelE; ntraj = $Ntraj; nstep = $Nstep; dt = 0.5
            tmax = 200; ntasks = $Ntasks; back_replicas = $backReplicas
            job_id = $jobId; server_dir = $runDir
        })
        Write-Output "$caseId v061 repeat=$repeat back=$backReplicas job=$jobId"
    }

    $qmCaseId = "qm_$($candidate.param_id)"
    if (-not (Already-Submitted $qmCaseId "qm" 0)) {
        $jobName = "qmp{0:d2}" -f $candidateNumber
        $remote = "cd $RemoteRoot && sbatch --parsable --ntasks=1 " +
            "--job-name=$jobName --export=ALL,AHM_WC_EV=$wc," +
            "AHM_ETA=$eta,AHM_DELE_EV=$DelE,AHM_NSTEP=$Nstep " +
            "run_mpi_cloud.slurm 100 10 5 0"
        $jobId = Invoke-Remote $remote
        Save-Row ([pscustomobject]@{
            stage = "pilot"; case_id = $qmCaseId; param_id = $candidate.param_id
            origin = $candidate.origin; random_seed = $candidate.random_seed
            method = "qm"; repeat = 0; wc_eV = $wc; eta = $eta
            delE_eV = $DelE; ntraj = 0; nstep = $Nstep; dt = 0.5
            tmax = 200; ntasks = 1; back_replicas = 0; job_id = $jobId
            server_dir = "$RemoteRoot/cloud-runs/$jobId"
        })
        Write-Output "$qmCaseId qm job=$jobId"
    }
}

Write-Output "Pilot candidates: $CandidateFile"
Write-Output "Pilot submission table: $SubmissionFile"
