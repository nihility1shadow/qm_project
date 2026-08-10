param(
    [Parameter(Mandatory = $true)]
    [ValidateSet("refine", "fine", "validate", "confirm", "confirm2")]
    [string]$Stage,

    [Parameter(Mandatory = $true)]
    [string]$CandidateFile,

    [ValidateRange(1, 24000000)]
    [int]$Ntraj = 2000000,

    [ValidateRange(2, 5)]
    [int]$Repeats = 3,

    [ValidateRange(1, 256)]
    [int]$Ntasks = 128
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = "E:\project_qm"
$ScanRoot = Join-Path $ProjectRoot "scan_v061_oxygen_t200_q10_20260810"
$CandidatePath = if ([System.IO.Path]::IsPathRooted($CandidateFile)) {
    $CandidateFile
} else {
    Join-Path $ProjectRoot $CandidateFile
}
$SubmissionFile = Join-Path $ScanRoot "${Stage}_submission.csv"
$RemoteRoot = "/data/home/yd101802/yd101802/nonadia"
$RemoteScanRoot = "$RemoteRoot/cloud-runs/v061-oxygen-t200-q10-20260810/$Stage"
$Key = Join-Path $env:TEMP "qm_codex_ssh\qm_codex_ed25519"
$KnownHosts = Join-Path $ProjectRoot "tmp\cloud_known_hosts_20260809"
$HostName = "yd101802@36.212.18.57"
$DelE = -17.194874968839816
$Nstep = 400

$candidates = @(Import-Csv $CandidatePath)
if ($candidates.Count -lt 1) { throw "No candidates in $CandidatePath" }
foreach ($candidate in $candidates) {
    $wc = [double]$candidate.wc_eV
    $eta = [double]$candidate.eta
    $backReplicas = [int]$candidate.back_replicas
    $stratifyForward = [int]$candidate.stratify_forward
    if ($wc -le 0.0 -or $wc -gt 10.0) { throw "wc outside (0, 10] eV: $wc" }
    if ($eta -le 0.0) { throw "eta must be positive: $eta" }
    if ($backReplicas -lt 1 -or $backReplicas -gt 64) {
        throw "back_replicas outside [1, 64]: $backReplicas"
    }
    if ($stratifyForward -notin @(0, 1)) {
        throw "stratify_forward must be 0 or 1: $stratifyForward"
    }
}

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
    $stratifyForward = [int]$candidate.stratify_forward

    foreach ($repeat in 1..$Repeats) {
        if (Already-Submitted $caseId "v061" $repeat) { continue }
        $runDir = "$RemoteScanRoot/$caseId/rep$repeat"
        $jobName = "v61$($Stage.Substring(0,1)){0:d2}r$repeat" -f $candidateNumber
        $wrap = "set -euo pipefail; module purge; module load gcc/12.1.0; " +
            "module load openmpi/4.1.8; export OMP_NUM_THREADS=1 " +
            "MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; " +
            "export LD_LIBRARY_PATH=$RemoteRoot/.deps/install/lib:" +
            '$LD_LIBRARY_PATH' + "; mpirun --bind-to core --map-by core " +
            "-np $Ntasks $RemoteRoot/na_mpi_v061.out 10 1 $Ntraj 10 5 " +
            "> program.out 2>&1"
        $remote = "cd $RemoteRoot && mkdir -p $runDir && " +
            "sbatch --parsable --partition=cmh --ntasks=$Ntasks " +
            "--cpus-per-task=1 --mem-per-cpu=1G --time=12:00:00 " +
            "--job-name=$jobName --chdir=$runDir " +
            "--export=ALL,AHM_WC_EV=$wc,AHM_ETA=$eta,AHM_DELE_EV=$DelE," +
            "AHM_NSTEP=$Nstep,SEP_MB_BACK_REPLICAS=$backReplicas," +
            "SEP_MB_STRATIFY_FORWARD_COUNT=$stratifyForward " +
            "--output=slurm-%j.out --error=slurm-%j.err --wrap='$wrap'"
        $jobId = Invoke-Remote $remote
        Save-Row ([pscustomobject]@{
            stage = $Stage; case_id = $caseId; param_id = $candidate.param_id
            origin = $candidate.origin; random_seed = $candidate.random_seed
            method = "v061"; repeat = $repeat; wc_eV = $wc; eta = $eta
            delE_eV = $DelE; ntraj = $Ntraj; nstep = $Nstep; dt = 0.5
            tmax = 200; ntasks = $Ntasks; back_replicas = $backReplicas
            stratify_forward = $stratifyForward; job_id = $jobId
            server_dir = $runDir
        })
        Write-Output "$caseId v061 repeat=$repeat job=$jobId"
    }

    if (-not (Already-Submitted $caseId "qm" 0)) {
        $jobName = "qm$($Stage.Substring(0,1)){0:d2}" -f $candidateNumber
        $remote = "cd $RemoteRoot && sbatch --parsable --ntasks=1 " +
            "--job-name=$jobName --export=ALL,AHM_WC_EV=$wc," +
            "AHM_ETA=$eta,AHM_DELE_EV=$DelE,AHM_NSTEP=$Nstep " +
            "run_mpi_cloud.slurm 100 10 5 0"
        $jobId = Invoke-Remote $remote
        Save-Row ([pscustomobject]@{
            stage = $Stage; case_id = $caseId; param_id = $candidate.param_id
            origin = $candidate.origin; random_seed = $candidate.random_seed
            method = "qm"; repeat = 0; wc_eV = $wc; eta = $eta
            delE_eV = $DelE; ntraj = 0; nstep = $Nstep; dt = 0.5
            tmax = 200; ntasks = 1; back_replicas = 0
            stratify_forward = 0; job_id = $jobId
            server_dir = "$RemoteRoot/cloud-runs/$jobId"
        })
        Write-Output "$caseId qm job=$jobId"
    }
}

Write-Output "Submitted stage=$Stage candidates=$($candidates.Count) ntraj=$Ntraj repeats=$Repeats"
