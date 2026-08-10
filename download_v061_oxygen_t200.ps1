param(
    [ValidateSet("pilot", "screen", "refine", "fine", "validate", "confirm", "confirm2", "confirm3", "confirm4")]
    [string]$Stage = "pilot"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = "E:\project_qm"
$ScanRoot = Join-Path $ProjectRoot "scan_v061_oxygen_t200_q10_20260810"
$SubmissionFile = Join-Path $ScanRoot "${Stage}_submission.csv"
$DataRoot = Join-Path $ScanRoot "${Stage}_data"
$RawRoot = Join-Path $ScanRoot "${Stage}_raw"
$Key = Join-Path $env:TEMP "qm_codex_ssh\qm_codex_ed25519"
$KnownHosts = Join-Path $ProjectRoot "tmp\cloud_known_hosts_20260809"
$HostName = "yd101802@36.212.18.57"

New-Item -ItemType Directory -Force -Path $DataRoot, $RawRoot | Out-Null
$rows = Import-Csv $SubmissionFile

function Copy-RemoteTree([string]$RemotePath, [string]$LocalPath) {
    for ($attempt = 1; $attempt -le 8; $attempt++) {
        $previousPreference = $ErrorActionPreference
        $ErrorActionPreference = "Continue"
        & scp -q -r -i $Key -o BatchMode=yes -o UserKnownHostsFile=$KnownHosts `
            -P 22 "${HostName}:$RemotePath" $LocalPath 2>&1 | Out-Null
        $scpExitCode = $LASTEXITCODE
        $ErrorActionPreference = $previousPreference
        if ($scpExitCode -eq 0) { return }
        if ($attempt -lt 8) { Start-Sleep -Seconds 3 }
    }
    throw "Download failed after 8 attempts: $RemotePath"
}

$poissonRaw = Join-Path $RawRoot $Stage
if (-not (Test-Path $poissonRaw)) {
    $remoteStageRoot = "/data/home/yd101802/yd101802/nonadia/cloud-runs/" +
        "v061-oxygen-t200-q10-20260810/$Stage"
    Copy-RemoteTree $remoteStageRoot $RawRoot
}

$qmRows = @($rows | Where-Object { $_.method -eq "qm" })
foreach ($row in $qmRows) {
    $localQmRoot = Join-Path $RawRoot "qm-$($row.job_id)"
    if (-not (Test-Path $localQmRoot)) {
        Copy-RemoteTree $row.server_dir $RawRoot
        Rename-Item -LiteralPath (Join-Path $RawRoot $row.job_id) `
            -NewName "qm-$($row.job_id)"
    }
}

$qmByParam = @{}
foreach ($row in $qmRows) { $qmByParam[$row.param_id] = $row }

foreach ($row in @($rows | Where-Object { $_.method -eq "v061" })) {
    $destination = Join-Path $DataRoot "$($row.case_id)\v061\rep$($row.repeat)"
    New-Item -ItemType Directory -Force -Path $destination | Out-Null
    $dataName = "ahm-sepmb-s10-n5-$($row.ntraj).dat"
    $source = Join-Path $poissonRaw "$($row.case_id)\rep$($row.repeat)"
    if (-not (Test-Path $source)) {
        $sourceParent = Split-Path $source -Parent
        New-Item -ItemType Directory -Force -Path $sourceParent | Out-Null
        Copy-RemoteTree $row.server_dir $sourceParent
    }
    foreach ($name in @($dataName, "program.out")) {
        $sourceFile = Join-Path $source $name
        if (-not (Test-Path $sourceFile)) { throw "Missing file: $sourceFile" }
        Copy-Item -LiteralPath $sourceFile -Destination $destination -Force
    }
}

foreach ($paramId in $qmByParam.Keys) {
    $row = $qmByParam[$paramId]
    $destination = Join-Path $DataRoot "qm\$paramId"
    New-Item -ItemType Directory -Force -Path $destination | Out-Null
    $source = Join-Path $RawRoot "qm-$($row.job_id)"
    foreach ($name in @("ahm-qm-s10-n5.dat", "program.out")) {
        $sourceFile = Join-Path $source $name
        if (-not (Test-Path $sourceFile)) { throw "Missing file: $sourceFile" }
        Copy-Item -LiteralPath $sourceFile -Destination $destination -Force
    }
}

Write-Output "Downloaded and curated data: $DataRoot"
