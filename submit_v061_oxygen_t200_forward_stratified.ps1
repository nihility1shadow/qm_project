$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = "E:\project_qm"
$ScanRoot = Join-Path $ProjectRoot "scan_v061_oxygen_t200_q10_20260810"
$SubmissionFile = Join-Path $ScanRoot "pilot_submission.csv"
$RemoteRoot = "/data/home/yd101802/yd101802/nonadia"
$RemoteScanRoot = "$RemoteRoot/cloud-runs/v061-oxygen-t200-q10-20260810/pilot"
$Key = Join-Path $env:TEMP "qm_codex_ssh\qm_codex_ed25519"
$KnownHosts = Join-Path $ProjectRoot "tmp\cloud_known_hosts_20260809"
$HostName = "yd101802@36.212.18.57"
$CaseId = "p08_random01_stratf"
$ParamId = "random01"
$Wc = 4.145914
$Eta = 0.000124380
$DelE = -17.194874968839816
$Ntraj = 8000000
$Nstep = 400
$Ntasks = 128
$BackReplicas = 12

$rows = [System.Collections.Generic.List[object]]::new()
foreach ($row in Import-Csv $SubmissionFile) { $rows.Add($row) }
foreach ($repeat in 1..3) {
    $duplicate = $rows | Where-Object {
        $_.case_id -eq $CaseId -and $_.method -eq "v061" -and
        [int]$_.repeat -eq $repeat
    }
    if ($duplicate) { throw "Refusing duplicate submission: $CaseId repeat=$repeat" }
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

foreach ($repeat in 1..3) {
    $runDir = "$RemoteScanRoot/$CaseId/rep$repeat"
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
        "--job-name=v61sfr$repeat --chdir=$runDir " +
        "--export=ALL,AHM_WC_EV=$Wc,AHM_ETA=$Eta,AHM_DELE_EV=$DelE," +
        "AHM_NSTEP=$Nstep,SEP_MB_BACK_REPLICAS=$BackReplicas," +
        "SEP_MB_STRATIFY_FORWARD_COUNT=1 " +
        "--output=slurm-%j.out --error=slurm-%j.err --wrap='$wrap'"
    $jobId = Invoke-Remote $remote
    $rows.Add([pscustomobject]@{
        stage = "pilot"; case_id = $CaseId; param_id = $ParamId
        origin = "forward_stratified"; random_seed = 20260810
        method = "v061"; repeat = $repeat; wc_eV = $Wc; eta = $Eta
        delE_eV = $DelE; ntraj = $Ntraj; nstep = $Nstep; dt = 0.5
        tmax = 200; ntasks = $Ntasks; back_replicas = $BackReplicas
        job_id = $jobId; server_dir = $runDir
    })
    $rows | Export-Csv -NoTypeInformation -Encoding utf8 $SubmissionFile
    Write-Output "$CaseId repeat=$repeat stratify_forward=1 job=$jobId"
}
