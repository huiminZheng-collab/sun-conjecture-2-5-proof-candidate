$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$manifest = Join-Path $root 'RELEASE-MANIFEST.sha256'
if (-not (Test-Path -LiteralPath $manifest)) { throw 'Manifest missing.' }
$failed = $false
foreach ($line in Get-Content -LiteralPath $manifest -Encoding UTF8) {
    if ($line -notmatch '^([0-9a-f]{64})  (.+)$') {
        throw "Malformed manifest line: $line"
    }
    $expected = $matches[1]
    $relative = $matches[2].Replace('/', [IO.Path]::DirectorySeparatorChar)
    $path = Join-Path $root $relative
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        Write-Error "Missing: $relative"
        $failed = $true
        continue
    }
    $actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant()
    if ($actual -ne $expected) {
        Write-Error "Hash mismatch: $relative"
        $failed = $true
    }
}
if ($failed) { exit 1 }
Write-Output 'PASS: every manifested file exists and matches its SHA-256 hash.'
