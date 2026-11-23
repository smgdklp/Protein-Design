# statgen 合成样本回归测试
# 用法: powershell -ExecutionPolicy Bypass -File run_tests.ps1 [-Exe <statgen.exe>] [-Keep]

param(
    [string]$Exe = "",
    [switch]$Keep
)

$ErrorActionPreference = 'Stop'
$inv = [System.Globalization.CultureInfo]::InvariantCulture

$here = Split-Path -Parent $MyInvocation.MyCommand.Path
if ([string]::IsNullOrEmpty($Exe)) {
    $Exe = Join-Path $here '..\..\..\bin\Release\statgen.exe'
}
if (-not (Test-Path -LiteralPath $Exe)) { throw "statgen.exe not found: $Exe" }
$Exe = (Resolve-Path -LiteralPath $Exe).Path

$repoRoot = (Resolve-Path -LiteralPath (Join-Path $here '..\..\..')).Path
$workRoot = Join-Path $repoRoot 'build\test_runs'
New-Item -ItemType Directory -Path $workRoot -Force | Out-Null
$workDir = Join-Path $workRoot ('stat_' + [guid]::NewGuid().ToString('N'))
$dataDir = Join-Path $workDir 'data'
$outDir = Join-Path $workDir 'out'
New-Item -ItemType Directory -Path $dataDir -Force | Out-Null
New-Item -ItemType Directory -Path $outDir -Force | Out-Null

function Fmt([double]$v) { return $v.ToString('F6', $inv) }

function SegmentLine([string]$seg, [string]$seq, [double[]]$cp, [double[]]$sp, [double[]]$cq, [double[]]$sq) {
    $cpj = ($cp | ForEach-Object { Fmt $_ }) -join ','
    $spj = ($sp | ForEach-Object { Fmt $_ }) -join ','
    $cqj = ($cq | ForEach-Object { Fmt $_ }) -join ','
    $sqj = ($sq | ForEach-Object { Fmt $_ }) -join ','
    $n = $seq.Length
    $json = '{"v":1,"seg":"' + $seg + '","pdb":"T","chain":"A","start":1,"end":' + $n +
            ',"len":' + $n + ',"seq":"' + $seq + '","cos_phi":[' + $cpj + '],"sin_phi":[' + $spj +
            '],"cos_psi":[' + $cqj + '],"sin_psi":[' + $sqj + '],"src":"T.pdb"}'
    return $json
}

# 角度常量
$c60 = 0.5; $s60 = -0.866025        # phi = -60
$c45 = 0.707107; $s45 = -0.707107   # psi = -45

$lines = @()

# 1) AAA, phi=-60 psi=-45 -> bin(24,27)
$lines += SegmentLine 'T1_A_1' 'AAA' @(0, $c60, 0) @(0, $s60, 0) @($c45, $c45, 0) @($s45, $s45, 0)

# 2) AAA, phi=180 psi=0 -> phi 回绕到 bin 0, psi bin 36
$lines += SegmentLine 'T2_A_1' 'AAA' @(0, -1, 0) @(0, 0, 0) @(1, 1, 0) @(0, 0, 0)

# 3) AAAA, 两个 i 都是 (-60,-45)
$lines += SegmentLine 'T3_A_1' 'AAAA' @(0, $c60, $c60, 0) @(0, $s60, $s60, 0) @($c45, $c45, $c45, 0) @($s45, $s45, $s45, 0)

# 4) ACD, phi=0 psi=0 -> bin(36,36)
$lines += SegmentLine 'T4_A_1' 'ACD' @(0, 1, 0) @(0, 0, 0) @(1, 1, 0) @(0, 0, 0)

# 5) ACD, phi=-180 psi=180 -> 两个都回绕到 bin 0
$lines += SegmentLine 'T5_A_1' 'ACD' @(0, -1, 0) @(0, -0.0, 0) @(-1, -1, 0) @(0, -0.0, 0)

Set-Content -LiteralPath (Join-Path $dataDir 'segments.jsonl') -Value $lines -Encoding ascii

$cfgPath = Join-Path $workDir 'config.json'
$cfg = @{ BIN_PATH = $dataDir; RESULT_PATH = $outDir } | ConvertTo-Json
Set-Content -LiteralPath $cfgPath -Value $cfg -Encoding ascii

$inspectFile = Join-Path $workDir 'aaa.csv'
& $Exe --config $cfgPath --inspect A A A > $inspectFile
if ($LASTEXITCODE -ne 0) { throw "statgen exited with $LASTEXITCODE" }

$script:failures = 0
function Assert-That([bool]$Condition, [string]$Message) {
    if ($Condition) { Write-Host "PASS  $Message" } else { Write-Host "FAIL  $Message" -ForegroundColor Red; $script:failures++ }
}

# --- summary ---
$summary = @{}
Get-Content -LiteralPath (Join-Path $outDir 'summary.txt') | ForEach-Object {
    $kv = $_ -split '=', 2
    if ($kv.Count -eq 2) { $summary[$kv[0]] = $kv[1] }
}
Assert-That ($summary['samples'] -eq '6') "samples=6 (got $($summary['samples']))"
Assert-That ($summary['expected_samples'] -eq '6') "expected_samples=6"
Assert-That ($summary['nonempty_matrices'] -eq '2') "2 non-empty matrices (AAA, ACD)"
Assert-That ($summary['empty_matrices'] -eq '7998') "7998 empty matrices"
Assert-That ($summary['result'] -eq 'OK') "result OK"

# --- AAA 矩阵 ---
$rows = @()
$text = Get-Content -LiteralPath $inspectFile
foreach ($l in $text[1..72]) { $rows += ,($l -split ',') }
Assert-That ($rows.Count -eq 72) 'AAA matrix has 72 rows'
$sum = 0.0
foreach ($r in $rows) { foreach ($v in $r) { $sum += [double]$v } }
Assert-That ([math]::Abs($sum - 1.0) -lt 1e-5) "AAA matrix sums to 1 (got $sum)"
Assert-That ([math]::Abs([double]$rows[24][27] - 0.75) -lt 1e-6) "AAA P(24,27)=0.75"
Assert-That ([math]::Abs([double]$rows[0][36] - 0.25) -lt 1e-6) "AAA P(0,36)=0.25 (phi=+180 wrap)"

# --- ACD 矩阵 ---
$acdFile = Join-Path $workDir 'acd.csv'
& $Exe --config $cfgPath --inspect A C D > $acdFile
$rows2 = @()
$text2 = Get-Content -LiteralPath $acdFile
foreach ($l in $text2[1..72]) { $rows2 += ,($l -split ',') }
Assert-That ([math]::Abs([double]$rows2[36][36] - 0.5) -lt 1e-6) "ACD P(36,36)=0.5"
Assert-That ([math]::Abs([double]$rows2[0][0] - 0.5) -lt 1e-6) "ACD P(0,0)=0.5 (phi=-180/psi=180 wrap)"

# --- npy 头 ---
$npy = Join-Path $outDir 'triplet_probs.npy'
$fs = [System.IO.File]::OpenRead($npy)
$buf = New-Object byte[] 128
$null = $fs.Read($buf, 0, 128)
$fs.Close()
$head = [System.Text.Encoding]::ASCII.GetString($buf)
Assert-That ($buf[0] -eq 0x93 -and $buf[1] -eq 0x4E -and $buf[2] -eq 0x55) 'npy magic correct'
Assert-That ($head -match "'descr': '<f4'") 'npy descr <f4'
Assert-That ($head -match "\(8000, 72, 72\)") 'npy shape (8000, 72, 72)'

if (-not $Keep) { Remove-Item -LiteralPath $workDir -Recurse -Force -ErrorAction SilentlyContinue }
else { Write-Host "kept: $workDir" }

if ($script:failures -gt 0) { Write-Host "$($script:failures) test(s) FAILED" -ForegroundColor Red; exit 1 }
Write-Host 'ALL TESTS PASSED' -ForegroundColor Green
exit 0
