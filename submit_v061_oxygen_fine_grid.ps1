$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = "E:\project_qm"
$OutputRoot = Join-Path $ProjectRoot "scan_v061_oxygen_grid_20260809"
$SubmissionFile = Join-Path $OutputRoot "fine_submission.csv"
$RemoteRoot = "/data/home/yd101802/yd101802/nonadia"
$RemoteGridRoot = "$RemoteRoot/cloud-runs/v061-oxygen-grid-20260809/fine"
$Key = Join-Path $env:TEMP "qm_codex_ssh\qm_codex_ed25519"
$KnownHosts = Join-Path $ProjectRoot "tmp\cloud_known_hosts_20260809"
$HostName = "yd101802@36.212.18.57"
$DelE = -17.194874968839816
$Ntraj = 1000000
$Nstep = 150
$Ntasks = 32
$Repeats = 3
$WcValues = @(0.25, 0.5, 1.0, 2.0, 4.0)
$EtaValues = @(2e-4, 3e-4, 5e-4, 7e-4, 1e-3)

$cases = [System.Collections.Generic.List[object]]::new()
foreach ($wc in $WcValues) {
    foreach ($eta in $EtaValues) {
        $cases.Add([pscustomobject]@{ wc = $wc; eta = $eta })
    }
}
$cases.Add([pscustomobject]@{ wc = 7.0; eta = 3e-4 })
$cases.Add([pscustomobject]@{ wc = 10.0; eta = 3e-4 })

New-Item -ItemType Directory -Force -Path $OutputRoot | Out-Null
$rows = [System.Collections.Generic.List[object]]::new()
if (Test-Path $SubmissionFile) {
    foreach ($row in Import-Csv $SubmissionFile) { $rows.Add($row) }
}

function Format-Token([double]$Value) {
    return $Value.ToString("0.###############", [Globalization.CultureInfo]::InvariantCulture).
        Replace("-", "m").Replace(".", "p")
}

function Invoke-Remote([string]$Command) {
    for ($attempt = 1; $attempt -le 90; $attempt++) {
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
        if ($attempt -lt 90) {
            Start-Sleep -Seconds 10
            continue
        }
        throw "Submission failed after $attempt attempts"
    }
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

$caseNumber = 0
foreach ($case in $cases) {
    $caseNumber++
    $wc = [double]$case.wc
    $eta = [double]$case.eta
    $caseId = "f{0:d2}_wc{1}_eta{2}" -f $caseNumber,
        (Format-Token $wc), (Format-Token $eta)

    foreach ($repeat in 1..$Repeats) {
        if (Already-Submitted $caseId "v061" $repeat) { continue }
        $runDir = "$RemoteGridRoot/$caseId/rep$repeat"
        $jobName = "v61f{0:d2}r$repeat" -f $caseNumber
        $wrap = "set -euo pipefail; module purge; module load gcc/12.1.0; " +
            "module load openmpi/4.1.8; export OMP_NUM_THREADS=1 " +
            "MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; " +
            "export LD_LIBRARY_PATH=$RemoteRoot/.deps/install/lib:" +
            '$LD_LIBRARY_PATH' + "; mpirun --bind-to core --map-by core " +
            "-np $Ntasks $RemoteRoot/na_mpi_v061.out 10 1 $Ntraj 10 5 " +
            "> program.out 2>&1"
        $remote = "cd $RemoteRoot && mkdir -p $runDir && " +
            "sbatch --parsable --partition=cmh --nodes=1 --ntasks=$Ntasks " +
            "--cpus-per-task=1 --mem-per-cpu=1G --time=02:00:00 " +
            "--job-name=$jobName --chdir=$runDir " +
            "--export=ALL,AHM_WC_EV=$wc,AHM_ETA=$eta,AHM_DELE_EV=$DelE," +
            "AHM_NSTEP=$Nstep,SEP_MB_BACK_REPLICAS=4," +
            "SEP_MB_STRATIFY_FORWARD_COUNT=0 " +
            "--output=slurm-%j.out --error=slurm-%j.err --wrap='$wrap'"
        $jobId = Invoke-Remote $remote
        Save-Row ([pscustomobject]@{
            stage = "fine"; case_id = $caseId; method = "v061"
            repeat = $repeat; wc_eV = $wc; eta = $eta; delE_eV = $DelE
            ntraj = $Ntraj; nstep = $Nstep; dt = 0.5; tmax = 75
            ntasks = $Ntasks; back_replicas = 4; job_id = $jobId
            server_dir = $runDir
        })
        Write-Output "$caseId v061 repeat=$repeat job=$jobId"
    }

    if (-not (Already-Submitted $caseId "qm" 0)) {
        $jobName = "qmf{0:d2}" -f $caseNumber
        $remote = "cd $RemoteRoot && sbatch --parsable --ntasks=1 " +
            "--job-name=$jobName --export=ALL,AHM_WC_EV=$wc," +
            "AHM_ETA=$eta,AHM_DELE_EV=$DelE,AHM_NSTEP=$Nstep " +
            "run_mpi_cloud.slurm 100 10 5 0"
        $jobId = Invoke-Remote $remote
        Save-Row ([pscustomobject]@{
            stage = "fine"; case_id = $caseId; method = "qm"
            repeat = 0; wc_eV = $wc; eta = $eta; delE_eV = $DelE
            ntraj = 0; nstep = $Nstep; dt = 0.5; tmax = 75
            ntasks = 1; back_replicas = 0; job_id = $jobId
            server_dir = "$RemoteRoot/cloud-runs/$jobId"
        })
        Write-Output "$caseId qm job=$jobId"
    }
}

Write-Output "Fine submission table: $SubmissionFile"
