param(
    [string]$Python = 'python',
    [string]$DwsimBin = '',
    [switch]$SkipDwsim
)
$ErrorActionPreference = 'Stop'
$packageRoot = $PSScriptRoot
if (Test-Path -LiteralPath $Python -PathType Leaf) {
    $pythonExe = (Resolve-Path -LiteralPath $Python).Path
} else {
    $pythonExe = (Get-Command $Python -ErrorAction Stop).Source
}
$previousDwsimBin = $env:DWSIM_BIN
if ($DwsimBin) { $env:DWSIM_BIN = $DwsimBin }
$logDirectory = Join-Path $packageRoot 'results/reproduction_logs'
New-Item -ItemType Directory -Path $logDirectory -Force | Out-Null
Push-Location $packageRoot
try {
    $scriptNames = @('run_analysis.py', 'run_process_analysis.py')
    if (-not $SkipDwsim) { $scriptNames += @('run_dwsim.py', 'run_connected_dwsim.py') }
    foreach ($scriptName in $scriptNames) {
        Write-Host "Running $scriptName"
        $scriptFile = Join-Path $packageRoot "code/$scriptName"
        $logFile = Join-Path $logDirectory ($scriptName + '.log')
        & $pythonExe $scriptFile *> $logFile
        if ($LASTEXITCODE -ne 0) {
            Get-Content -LiteralPath $logFile -Tail 35
            throw "$scriptName failed; full output: $logFile"
        }
    }
    Write-Host 'Reproduction completed. Results and native verification reports are in the results folder.'
} finally {
    Pop-Location
    if ($null -eq $previousDwsimBin) {
        Remove-Item Env:DWSIM_BIN -ErrorAction SilentlyContinue
    } else {
        $env:DWSIM_BIN = $previousDwsimBin
    }
}
