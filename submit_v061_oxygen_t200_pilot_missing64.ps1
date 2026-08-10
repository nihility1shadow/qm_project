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

$missing = @(
    [pscustomobject]@{ case_id = "p02_center_b08"; repeat = 3 },
    [pscustomobject]@{ case_id = "p03_center_b12"; repeat = 1 },
    [pscustomobject]@{ case_id = "p03_center_b12"; repeat = 2 },
    [pscustomobject]@{ case_id = "p03_center_b12"; repeat = 3 },
    [pscustomobject]@{ case_id = "p04_center_b16"; repeat = 1 },
    [pscustomobject]@{ case_id = "p04_center_b16"; repeat = 2 },
    [pscustomobject]@{ case_id = "p04_center_b16"; repeat = 3 },
    [pscustomobject]@{ case_id = "p05_random_01"; repeat = 1 },
    [pscustomobject]@{ case_id = "p05_random_01"; repeat = 2 },
    [pscustomobject]@{ case_id = "p05_random_01"; repeat = 3 },
    [pscustomobject]@{ case_id = "p06_random_02"; repeat = 1 },
    [pscustomobject]@{ case_id = "p06_random_02"; repeat = 2 },
    [pscustomobject]@{ case_id = "p06_random_02"; repeat = 3 },
    [pscustomobject]@{ case_id = "p07_random_03"; repeat = 1 },
    [pscustomobject]@{ case_id = "p07_random_03"; repeat = 2 },
    [pscustomobject]@{ case_id = "p07_random_03"; repeat = 3 }
)

$candidates = @{}
foreach ($candidate in Import-Csv $CandidateFile) {
    $candidates[$candidate.case_id] = $candidate
}
$rows = [System.Collections.Generic.List[object]]::new()
foreach ($row in Import-Csv $SubmissionFile) { $rows.Add($row) }

foreach ($item in $missing) {
    $duplicate = $rows | Where-Object {
        $_.case_id -eq $item.case_id -and $_.method -eq "v061" -and
        [int]$_.repeat -eq [int]$item.repeat
    }
    if ($duplicate) {
        throw "Refusing duplicate submission: $($item.case_id) repeat=$($item.repeat)"
    }
    if (-not $candidates.ContainsKey($item.case_id)) {
        throw "Unknown candidate: $($item.case_id)"
    }
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

foreach ($item in $missing) {
    $candidate = $candidates[$item.case_id]
    $caseId = $candidate.case_id
    $repeat = [int]$item.repeat
    $wc = [double]$candidate.wc_eV
    $eta = [double]$candidate.eta
    $backReplicas = [int]$candidate.back_replicas
    $runDir = "$RemoteScanRoot/$caseId/rep$repeat"
    $jobName = "v61m$repeat"
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
    $row = [pscustomobject]@{
        stage = "pilot"; case_id = $caseId; param_id = $candidate.param_id
        origin = $candidate.origin; random_seed = $candidate.random_seed
        method = "v061"; repeat = $repeat; wc_eV = $wc; eta = $eta
        delE_eV = $DelE; ntraj = $Ntraj; nstep = $Nstep; dt = 0.5
        tmax = 200; ntasks = $Ntasks; back_replicas = $backReplicas
        job_id = $jobId; server_dir = $runDir
    }
    $rows.Add($row)
    $rows | Export-Csv -NoTypeInformation -Encoding utf8 $SubmissionFile
    Write-Output "$caseId repeat=$repeat back=$backReplicas job=$jobId"
}

Write-Output "Submitted exactly $($missing.Count) whitelisted jobs."
