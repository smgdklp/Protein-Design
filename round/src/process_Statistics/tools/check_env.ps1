# uv 工具链检测（Windows PowerShell 5.1 兼容）
# 用法: powershell -ExecutionPolicy Bypass -File check_env.ps1 [-Install]
#   -Install  在项目 .venv 里离线安装 numpy / matplotlib（走 uv 缓存，不联网）

param([switch]$Install)

$ErrorActionPreference = 'Stop'

$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$proj = Split-Path -Parent $here
$venv = Join-Path $proj '.venv'
$venvPy = Join-Path $venv 'Scripts\python.exe'

$uv = $null
if ($env:UV -and (Test-Path -LiteralPath $env:UV)) {
    $uv = $env:UV
} else {
    $candidates = @('C:\Users\ThinkPad\.local\bin\uv.exe')
    $cmd = Get-Command uv -ErrorAction SilentlyContinue
    if ($cmd) { $candidates += $cmd.Source }
    foreach ($c in $candidates) {
        if ($c -and (Test-Path -LiteralPath $c)) { $uv = $c; break }
    }
}
if (-not $uv) {
    Write-Host 'FAIL: uv not found (expected C:\Users\ThinkPad\.local\bin\uv.exe)' -ForegroundColor Red
    exit 1
}

$uvVersion = (& $uv --version 2>&1) -join ' '
Write-Host "uv          : $uv"
Write-Host "uv version  : $uvVersion"

$py312 = (& $uv python find 3.12 2>$null)
if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($py312)) {
    Write-Host 'FAIL: uv-managed CPython 3.12 not found' -ForegroundColor Red
    exit 1
}
Write-Host "python 3.12 : $py312"

if (-not (Test-Path -LiteralPath $venvPy)) {
    Write-Host "creating venv: $venv"
    & $uv venv --python 3.12 $venv
    if ($LASTEXITCODE -ne 0) { Write-Host 'FAIL: uv venv failed' -ForegroundColor Red; exit 1 }
}

if ($Install) {
    Write-Host 'installing numpy / matplotlib from uv cache (offline)'
    & $uv pip install --offline --python $venvPy numpy matplotlib
    if ($LASTEXITCODE -ne 0) {
        Write-Host 'FAIL: uv pip install --offline failed' -ForegroundColor Red
        exit 1
    }
}

$probe = @'
import sys, numpy, matplotlib
print(sys.prefix)
print(numpy.__version__)
print(matplotlib.__version__)
print(matplotlib.get_backend())
'@
$probeFile = Join-Path $env:TEMP ('stat_probe_' + [guid]::NewGuid().ToString('N') + '.py')
Set-Content -LiteralPath $probeFile -Value $probe -Encoding ascii
$out = & $venvPy $probeFile 2>&1
$code = $LASTEXITCODE
Remove-Item -LiteralPath $probeFile -Force -ErrorAction SilentlyContinue

if ($code -ne 0) {
    Write-Host 'FAIL: cannot import numpy / matplotlib in venv' -ForegroundColor Red
    $out | ForEach-Object { Write-Host "  $_" }
    exit 1
}

$linesOut = @($out)
Write-Host "prefix      : $($linesOut[0])"
Write-Host "numpy       : $($linesOut[1])"
Write-Host "matplotlib  : $($linesOut[2])"
Write-Host "backend     : $($linesOut[3])"

$cfgFile = Join-Path $venv 'pyvenv.cfg'
$uvManaged = $false
if (Test-Path -LiteralPath $cfgFile) {
    $uvManaged = ((Get-Content -LiteralPath $cfgFile) -match '^uv\s*=').Count -gt 0
}
Write-Host ("uv-managed  : " + $(if ($uvManaged) { 'yes' } else { 'NO' }))

if (-not $uvManaged) { Write-Host 'FAIL: venv is not uv-managed' -ForegroundColor Red; exit 1 }
if (-not ($linesOut[0] -like "$proj*")) { Write-Host 'FAIL: python prefix is not the project venv' -ForegroundColor Red; exit 1 }

Write-Host 'ENV OK' -ForegroundColor Green
exit 0
