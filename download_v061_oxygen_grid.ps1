param(
    [ValidateSet("coarse", "fine")]
    [string]$Stage = "coarse"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = "E:\project_qm"
$GridRoot = Join-Path $ProjectRoot "scan_v061_oxygen_grid_20260809"
$SubmissionFile = Join-Path $GridRoot "${Stage}_submission.csv"
$DataRoot = Join-Path $GridRoot "${Stage}_data"
$Key = Join-Path $env:TEMP "qm_codex_ssh\qm_codex_ed25519"
$KnownHosts = Join-Path $ProjectRoot "tmp\cloud_known_hosts_20260809"
$HostName = "yd101802@36.212.18.57"

New-Item -ItemType Directory -Force -Path $DataRoot | Out-Null
$rows = Import-Csv $SubmissionFile

function Copy-RemoteFile([string]$RemotePath, [string]$LocalPath) {
    for ($attempt = 1; $attempt -le 8; $attempt++) {
        $previousPreference = $ErrorActionPreference
        $ErrorActionPreference = "Continue"
        & scp -q -i $Key -o BatchMode=yes -o UserKnownHostsFile=$KnownHosts `
            -P 22 "${HostName}:$RemotePath" $LocalPath 2>&1 | Out-Null
        $scpExitCode = $LASTEXITCODE
        $ErrorActionPreference = $previousPreference
        if ($scpExitCode -eq 0 -and (Test-Path $LocalPath) -and
            (Get-Item $LocalPath).Length -gt 0) {
            return
        }
        Remove-Item -LiteralPath $LocalPath -Force -ErrorAction SilentlyContinue
        if ($attempt -lt 8) { Start-Sleep -Seconds 3 }
    }
    throw "Failed to download $RemotePath after 8 attempts"
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
    $localData = Join-Path $destination $dataName
    if (-not (Test-Path $localData)) {
        Copy-RemoteFile "$($row.server_dir)/$dataName" $localData
    }

    $localProgram = Join-Path $destination "program.out"
    if (-not (Test-Path $localProgram)) {
        Copy-RemoteFile "$($row.server_dir)/program.out" $localProgram
    }
    Write-Output "$($row.case_id) $($row.method) $repeatToken"
}

Write-Output "Downloaded data: $DataRoot"
