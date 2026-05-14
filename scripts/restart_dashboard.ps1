# Restart CASCADE FastAPI (python -m dashboard_backend applies Proactor on Windows).
# Usage: .\scripts\restart_dashboard.ps1   or   .\scripts\restart_dashboard.ps1 -Port 8000

param(
    [int]$Port = 8001,
    [switch]$Reload
)

$ErrorActionPreference = "Stop"
$root = Resolve-Path (Join-Path $PSScriptRoot "..")
Set-Location $root

foreach ($ln in @(Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue)) {
    Write-Host "Stopping PID $($ln.OwningProcess) on port $Port"
    Stop-Process -Id $ln.OwningProcess -Force -ErrorAction SilentlyContinue
}
Start-Sleep -Seconds 2

$pythonExe = Join-Path $root ".venv\Scripts\python.exe"
if (!(Test-Path $pythonExe)) {
    Write-Error "Missing $pythonExe -- create .venv first."
    exit 1
}

if (Test-Path (Join-Path $root ".env")) {
    Get-Content (Join-Path $root ".env") | ForEach-Object {
        if ($_ -match '^\s*([^#=]+)=(.*)$') {
            $k = $Matches[1].Trim()
            $v = $Matches[2].Trim().Trim('"').Trim("'")
            [System.Environment]::SetEnvironmentVariable($k, $v, "Process")
        }
    }
}

$userScripts = (& $pythonExe -c "import sysconfig,sys; sys.stdout.write(sysconfig.get_path('scripts','nt_user'))").Trim()
if (Test-Path (Join-Path $userScripts "vastai.exe")) {
    $env:PATH = "$userScripts;$env:PATH"
}

$argv = @("-m", "dashboard_backend", "--host", "127.0.0.1", "--port", "$Port")
if (-not $Reload) {
    $argv += "--no-reload"
}

$outLog = Join-Path $root "outputs\backend.log"
$errLog = Join-Path $root "outputs\backend.err.log"
New-Item -ItemType Directory -Path (Join-Path $root "outputs") -Force | Out-Null

$proc = Start-Process -FilePath $pythonExe -ArgumentList $argv -WorkingDirectory $root `
    -RedirectStandardOutput $outLog -RedirectStandardError $errLog -WindowStyle Hidden -PassThru

Write-Host "Started backend PID $($proc.Id) on http://127.0.0.1:$Port (logs outputs\backend*.log)"
