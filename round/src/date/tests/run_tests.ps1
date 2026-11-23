# PDB 片段提取器 —— 合成样本回归测试
# 用法: pwsh -File run_tests.ps1 [-Exe <pdbextract.exe>] [-Keep]

param(
    [string]$Exe = "",
    [switch]$Keep
)

$ErrorActionPreference = 'Stop'

$here = Split-Path -Parent $MyInvocation.MyCommand.Path
if ([string]::IsNullOrEmpty($Exe)) {
    $Exe = Join-Path $here '..\..\..\bin\Release\pdbextract.exe'
}
if (-not (Test-Path -LiteralPath $Exe)) {
    throw "pdbextract.exe not found: $Exe"
}
$Exe = (Resolve-Path -LiteralPath $Exe).Path

$repoRoot = (Resolve-Path -LiteralPath (Join-Path $here '..\..\..')).Path
$workRoot = Join-Path $repoRoot 'build\test_runs'
New-Item -ItemType Directory -Path $workRoot -Force | Out-Null
$workDir = Join-Path $workRoot ([guid]::NewGuid().ToString('N'))
$dataDir = Join-Path $workDir 'data'
$outDir = Join-Path $workDir 'out'
New-Item -ItemType Directory -Path $dataDir -Force | Out-Null
New-Item -ItemType Directory -Path $outDir -Force | Out-Null

$script:lines = $null
$script:serial = 0

function New-AtomLine {
    param(
        [string]$Atom,
        [string]$Res,
        [string]$Chain,
        [int]$Seq,
        [double]$X,
        [double]$Y,
        [double]$Z,
        [string]$Record = 'ATOM'
    )
    $script:serial++
    $rec = if ($Record -eq 'HETATM') { 'HETATM' } else { 'ATOM  ' }
    $s = $rec
    $s += ('{0,5}' -f $script:serial)
    $s += ' '
    $s += ('{0,-4}' -f $Atom)
    $s += ' '
    $s += ('{0,3}' -f $Res)
    $s += ' '
    $s += $Chain
    $s += ('{0,4}' -f $Seq)
    $s += ' '
    $s += '   '
    $s += ('{0,8:0.000}' -f $X)
    $s += ('{0,8:0.000}' -f $Y)
    $s += ('{0,8:0.000}' -f $Z)
    $s += ('{0,6:0.00}' -f 1.0)
    $s += ('{0,6:0.00}' -f 0.0)
    $s += '          '
    $s += ('{0,2}' -f 'C')
    return $s
}

function Add-Peptide {
    param(
        [string]$Res,
        [string]$Chain,
        [int]$StartSeq,
        [int]$Count,
        [double]$XOffset = 0.0,
        [switch]$SkipCA
    )
    for ($i = 0; $i -lt $Count; $i++) {
        $seq = $StartSeq + $i
        $x = $XOffset + $i * 3.8
        $script:lines.Add((New-AtomLine 'N' $Res $Chain $seq ($x - 1.2) 0.5 0.3))
        if (-not $SkipCA) {
            $script:lines.Add((New-AtomLine 'CA' $Res $Chain $seq $x 0.0 0.0))
        }
        $script:lines.Add((New-AtomLine 'C' $Res $Chain $seq ($x + 1.2) -0.4 0.7))
    }
}

function Write-TestPdb {
    param([string]$Name, [scriptblock]$Body)
    $script:lines = New-Object System.Collections.Generic.List[string]
    $script:serial = 0
    & $Body
    $script:lines.Add('END')
    $path = Join-Path $dataDir $Name
    Set-Content -LiteralPath $path -Value $script:lines -Encoding ascii
}

# --- 1. 纯标准 5 残基 ---
Write-TestPdb 'ok5.pdb' {
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 1 -Count 5
    $script:lines.Add('TER')
}

# --- 2. 只有 2 残基，应被最小长度 3 丢掉 ---
Write-TestPdb 'short2.pdb' {
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 1 -Count 2
    $script:lines.Add('TER')
}

# --- 3. 标准-非标准(DA)-标准，应切成两段 ---
Write-TestPdb 'mixed.pdb' {
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 1 -Count 4
    $script:lines.Add((New-AtomLine 'N' 'DA' 'A' 5 14.0 0.5 0.3))
    $script:lines.Add((New-AtomLine 'CA' 'DA' 'A' 5 15.2 0.0 0.0))
    $script:lines.Add((New-AtomLine 'C' 'DA' 'A' 5 16.4 -0.4 0.7))
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 6 -Count 4 -XOffset 19.0
    $script:lines.Add('TER')
}

# --- 4. MSE 保留并记为 MET ---
Write-TestPdb 'mse.pdb' {
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 1 -Count 1
    $script:lines.Add((New-AtomLine 'N' 'MSE' 'A' 2 2.6 0.5 0.3 -Record 'HETATM'))
    $script:lines.Add((New-AtomLine 'CA' 'MSE' 'A' 2 3.8 0.0 0.0 -Record 'HETATM'))
    $script:lines.Add((New-AtomLine 'C' 'MSE' 'A' 2 5.0 -0.4 0.7 -Record 'HETATM'))
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 3 -Count 1 -XOffset 7.6
    $script:lines.Add('TER')
}

# --- 5. 同链但坐标断裂，应切成两段 ---
Write-TestPdb 'far.pdb' {
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 1 -Count 3
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 4 -Count 3 -XOffset 20.0
    $script:lines.Add('TER')
}

# --- 6. 中间残基缺 CA，应切成两段 ---
Write-TestPdb 'missing.pdb' {
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 1 -Count 3
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 4 -Count 1 -XOffset 11.4 -SkipCA
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 5 -Count 3 -XOffset 15.2
    $script:lines.Add('TER')
}

# --- 7. 多模型，只取第一个 ---
Write-TestPdb 'multimodel.pdb' {
    $script:lines.Add('MODEL        1')
    Add-Peptide -Res 'ALA' -Chain 'A' -StartSeq 1 -Count 3
    $script:lines.Add('ENDMDL')
    $script:lines.Add('MODEL        2')
    Add-Peptide -Res 'GLY' -Chain 'A' -StartSeq 1 -Count 5
    $script:lines.Add('ENDMDL')
}

# --- 配置 ---
$cfgPath = Join-Path $workDir 'pbd_config.json'
$cfg = @{
    PBD_PATH = $dataDir
    BIN_PATH = $outDir
    NUM      = 2
} | ConvertTo-Json
Set-Content -LiteralPath $cfgPath -Value $cfg -Encoding ascii

& $Exe --config $cfgPath
if ($LASTEXITCODE -ne 0) {
    throw "pdbextract exited with $LASTEXITCODE"
}

# --- 断言 ---
$script:failures = 0
function Assert-That {
    param([bool]$Condition, [string]$Message)
    if ($Condition) {
        Write-Host "PASS  $Message"
    } else {
        Write-Host "FAIL  $Message" -ForegroundColor Red
        $script:failures++
    }
}

$segPath = Join-Path $outDir 'segments.jsonl'
$procPath = Join-Path $outDir 'processed.jsonl'
$segments = @()
if (Test-Path -LiteralPath $segPath) {
    $segments = @(Get-Content -LiteralPath $segPath | Where-Object { $_.Trim() } | ForEach-Object { $_ | ConvertFrom-Json })
}
$processed = @(Get-Content -LiteralPath $procPath | Where-Object { $_.Trim() } | ForEach-Object { $_ | ConvertFrom-Json })

function Segs-Of {
    param([string]$Id)
    return @($segments | Where-Object { $_.pdb -eq $Id })
}

$ok5 = @(Segs-Of 'OK5')
Assert-That ($ok5.Count -eq 1) 'ok5: 1 segment'
Assert-That ($ok5.Count -ge 1 -and $ok5[0].seq -eq 'AAAAA') 'ok5: sequence AAAAA'
Assert-That ($ok5.Count -ge 1 -and $ok5[0].len -eq 5) 'ok5: len 5'
Assert-That ($ok5.Count -ge 1 -and $ok5[0].chain -eq 'A') 'ok5: chain A'

Assert-That (@(Segs-Of 'SHORT2').Count -eq 0) 'short2: dropped by MIN_LEN=3'

$mixed = @(Segs-Of 'MIXED')
Assert-That ($mixed.Count -eq 2) 'mixed: 2 segments'
if ($mixed.Count -eq 2) {
    Assert-That ($mixed[0].seq -eq 'AAAA' -and $mixed[1].seq -eq 'AAAA') 'mixed: both AAAA'
}

$mse = @(Segs-Of 'MSE')
Assert-That ($mse.Count -eq 1) 'mse: 1 segment'
Assert-That ($mse.Count -eq 1 -and $mse[0].seq -eq 'AMA') 'mse: MSE mapped to M'

$far = @(Segs-Of 'FAR')
Assert-That ($far.Count -eq 2) 'far: 2 segments on geometry break'
if ($far.Count -eq 2) {
    Assert-That ($far[0].start -eq 1 -and $far[0].end -eq 3) 'far: first segment residues 1-3'
    Assert-That ($far[1].start -eq 4 -and $far[1].end -eq 6) 'far: second segment residues 4-6'
}

$missingSegs = @(Segs-Of 'MISSING')
Assert-That ($missingSegs.Count -eq 2) 'missing: 2 segments around missing CA'

$mm = @(Segs-Of 'MULTIMODEL')
Assert-That ($mm.Count -eq 1) 'multimodel: only first model'
Assert-That ($mm.Count -eq 1 -and $mm[0].seq -eq 'AAA') 'multimodel: sequence AAA'

Assert-That ($processed.Count -eq 7) "processed: 7 entries (got $($processed.Count))"
Assert-That (@($processed | Where-Object { $_.status -ne 'ok' }).Count -eq 0) 'processed: all ok'
Assert-That ($segments.Count -eq 9) "segments total 9 (got $($segments.Count))"

# 角度自洽：占位零点 + 其余 cos^2+sin^2=1
$angleErrors = 0
$placeholderErrors = 0
foreach ($s in $segments) {
    $n = $s.len
    if ($s.cos_phi.Count -ne $n -or $s.sin_phi.Count -ne $n -or
        $s.cos_psi.Count -ne $n -or $s.sin_psi.Count -ne $n) {
        $angleErrors++
        continue
    }
    if ([math]::Abs($s.cos_phi[0]) -gt 1e-9 -or [math]::Abs($s.sin_phi[0]) -gt 1e-9) { $placeholderErrors++ }
    if ([math]::Abs($s.cos_psi[$n - 1]) -gt 1e-9 -or [math]::Abs($s.sin_psi[$n - 1]) -gt 1e-9) { $placeholderErrors++ }
    for ($i = 0; $i -lt $n; $i++) {
        if ($i -ge 1) {
            $r = $s.cos_phi[$i] * $s.cos_phi[$i] + $s.sin_phi[$i] * $s.sin_phi[$i]
            if ([math]::Abs($r - 1.0) -gt 1e-4) { $angleErrors++ }
        }
        if ($i -le $n - 2) {
            $r = $s.cos_psi[$i] * $s.cos_psi[$i] + $s.sin_psi[$i] * $s.sin_psi[$i]
            if ([math]::Abs($r - 1.0) -gt 1e-4) { $angleErrors++ }
        }
    }
}
Assert-That ($angleErrors -eq 0) "angles: cos^2+sin^2==1 everywhere valid (errors=$angleErrors)"
Assert-That ($placeholderErrors -eq 0) "angles: terminal placeholders are zero (errors=$placeholderErrors)"

# 宽松 JSON：单反斜杠 + 尾逗号 + 字符串线程数
$sloppyDir = Join-Path $workDir 'sloppy'
New-Item -ItemType Directory -Path $sloppyDir -Force | Out-Null
$emptyIn = Join-Path $workDir 'empty_in'
New-Item -ItemType Directory -Path $emptyIn -Force | Out-Null
$sloppyCfg = Join-Path $workDir 'sloppy_config.json'
$sloppyText = @"
{
    "PBD_PATH" : "$emptyIn",
    "BIN_PATH" : "$sloppyDir",
    "NUM" : "2",
}
"@
Set-Content -LiteralPath $sloppyCfg -Value $sloppyText -Encoding ascii
& $Exe --config $sloppyCfg | Out-Null
Assert-That ($LASTEXITCODE -eq 0) 'sloppy config: parsed and exited 0'

# 统计模式
& $Exe --config $cfgPath --stats | Out-Null
Assert-That ($LASTEXITCODE -eq 0) '--stats runs on generated segments'

if (-not $Keep) {
    Remove-Item -LiteralPath $workDir -Recurse -Force -ErrorAction SilentlyContinue
} else {
    Write-Host "kept: $workDir"
}

if ($script:failures -gt 0) {
    Write-Host "$($script:failures) test(s) FAILED" -ForegroundColor Red
    exit 1
}
Write-Host 'ALL TESTS PASSED' -ForegroundColor Green
exit 0
