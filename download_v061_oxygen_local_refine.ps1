param(
    [ValidateSet("screen", "confirm")]
    [string]$Stage = "screen"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = "E:\project_qm"
$ScanRoot = Join-Path $ProjectRoot "scan_v061_oxygen_local_refine_20260810"
$SubmissionFile = Join-Path $ScanRoot "${Stage}_submission.csv"
$DataRoot = Join-Path $ScanRoot "${Stage}_data"
$RawRoot = Join-Path $ScanRoot "${Stage}_raw"
$Key = Join-Path $env:TEMP "qm_codex_ssh\qm_codex_ed25519"
$KnownHosts = Join-Path $ProjectRoot "tmp\cloud_known_hosts_20260809"
$HostName = "yd101802@36.212.18.57"

New-Item -ItemType Directory -Force -Path $DataRoot, $RawRoot | Out-Null
$rows = Import-Csv $SubmissionFile

function Copy-RemoteTree([string[]]$RemotePaths, [string]$LocalPath) {
    for ($attempt = 1; $attempt -le 8; $attempt++) {
        $previousPreference = $ErrorActionPreference
        $ErrorActionPreference = "Continue"
        $sources = $RemotePaths | ForEach-Object { "${HostName}:$_" }
        & scp -q -r -i $Key -o BatchMode=yes -o UserKnownHostsFile=$KnownHosts `
            -P 22 @sources $LocalPath 2>&1 | Out-Null
        $scpExitCode = $LASTEXITCODE
        $ErrorActionPreference = $previousPreference
        if ($scpExitCode -eq 0) { return }
        if ($attempt -lt 8) { Start-Sleep -Seconds 3 }
    }
    throw "Bulk download failed after 8 attempts"
}

$poissonRaw = Join-Path $RawRoot $Stage
if (-not (Test-Path $poissonRaw)) {
    $poissonParent = Split-Path $poissonRaw -Parent
    $remoteStageRoot = "/data/home/yd101802/yd101802/nonadia/cloud-runs/" +
        "v061-oxygen-local-refine-20260810/$Stage"
    Copy-RemoteTree @($remoteStageRoot) $poissonParent
}

$qmRaw = Join-Path $RawRoot "qm_jobs"
New-Item -ItemType Directory -Force -Path $qmRaw | Out-Null
$qmRows = @($rows | Where-Object { $_.method -eq "qm" })
$missingQm = @($qmRows | Where-Object {
    -not (Test-Path (Join-Path $qmRaw $_.job_id))
})
if ($missingQm.Count -gt 0) {
    Copy-RemoteTree @($missingQm | ForEach-Object { $_.server_dir }) $qmRaw
}

foreach ($row in $rows) {
    $repeatToken = if ($row.method -eq "qm") { "exact" } else { "rep$($row.repeat)" }
    $destination = Join-Path $DataRoot "$($row.case_id)\$($row.method)\$repeatToken"
    New-Item -ItemType Directory -Force -Path $destination | Out-Null
    $dataName = if ($row.method -eq "qm") {
        "ahm-qm-s10-n5.dat"
    } else {
        "ahm-sepmb-s10-n5-$($row.ntraj).dat"
    }

    $source = if ($row.method -eq "qm") {
        Join-Path $qmRaw $row.job_id
    } else {
        Join-Path $poissonRaw "$($row.case_id)\rep$($row.repeat)"
    }
    $localData = Join-Path $destination $dataName
    $localProgram = Join-Path $destination "program.out"
    if (-not (Test-Path (Join-Path $source $dataName))) {
        throw "Missing downloaded data: $source\$dataName"
    }
    Copy-Item -LiteralPath (Join-Path $source $dataName) -Destination $localData -Force
    Copy-Item -LiteralPath (Join-Path $source "program.out") -Destination $localProgram -Force
    Write-Output "$($row.case_id) $($row.method) $repeatToken"
}

Write-Output "Downloaded data: $DataRoot"
