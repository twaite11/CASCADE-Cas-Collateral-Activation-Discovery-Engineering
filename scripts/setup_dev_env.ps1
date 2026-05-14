# scripts/setup_dev_env.ps1
#
# One-shot local-dev bootstrap for the CASCADE dashboard on Windows.
#
# What it does (idempotent -- safe to re-run):
#   1. Installs the `vastai` Python CLI if missing.
#   2. Asks for your Vast.ai API key once, stores it via `vastai set api-key`
#      (writes ~/.config/vastai/vast_api_key) AND exports VAST_API_KEY into
#      the current shell.
#   3. Sets CASCADE_API_KEY=disabled for localhost dev (the backend will log
#      a one-time warning so you know auth is off).
#   4. Detects your default SSH key and exports VAST_SSH_KEY_PATH; if no key
#      exists it offers to create one.
#   5. Smoke-tests `vastai search offers` so you can confirm credentials work
#      before clicking around in the dashboard.
#
# Usage (from repo root, in PowerShell):
#   .\scripts\setup_dev_env.ps1
#
# After it finishes you must ALSO tell the browser UI the same dashboard key:
#   - Open the dashboard in your browser
#   - DevTools console:
#       localStorage.setItem("cascade_api_key", "disabled")
#   - Reload the page.

[CmdletBinding()]
param(
    [string]$VastApiKey,
    [string]$CascadeApiKey = "disabled",
    [string]$SshKeyPath
)

$ErrorActionPreference = "Stop"

function Write-Step($msg) { Write-Host "==> $msg" -ForegroundColor Cyan }
function Write-OK($msg)   { Write-Host "    $msg" -ForegroundColor Green }
function Write-Warn($msg) { Write-Host "    $msg" -ForegroundColor Yellow }
function Write-Err($msg)  { Write-Host "    $msg" -ForegroundColor Red }

function Get-PythonUserScriptsDir {
    # Use sysconfig.get_path('scripts', 'nt_user') -- gives the actual
    # %APPDATA%\Python\PythonNNN\Scripts dir.  site.USER_BASE returns the
    # parent on Python >=3.11, which would resolve to the wrong path.
    try {
        $dir = (& python -c "import sysconfig,sys; sys.stdout.write(sysconfig.get_path('scripts','nt_user'))" 2>$null)
        if ($LASTEXITCODE -ne 0 -or -not $dir) { return $null }
        return $dir.Trim()
    } catch {
        return $null
    }
}

# 1. vastai CLI
Write-Step "Checking vastai CLI..."
$vastCmd = Get-Command vastai -ErrorAction SilentlyContinue
if (-not $vastCmd) {
    Write-Warn "not found on PATH; installing via pip..."
    python -m pip install --upgrade vastai
    # pip installs into the per-user Scripts dir which is NOT on PATH by default
    # on Windows. Prepend it for this session so the rest of the script works.
    $userScripts = Get-PythonUserScriptsDir
    if ($userScripts -and (Test-Path (Join-Path $userScripts 'vastai.exe'))) {
        $env:PATH = "$userScripts;$env:PATH"
        Write-OK "Prepended $userScripts to PATH for this session"
        Write-Warn "To make it permanent, run this ONCE in a fresh PowerShell window:"
        Write-Host "        setx PATH `"$userScripts;`$env:PATH`"" -ForegroundColor DarkGray
        Write-Warn "(Or add the dir via System Properties > Environment Variables > PATH.)"
    }
    $vastCmd = Get-Command vastai -ErrorAction SilentlyContinue
    if (-not $vastCmd) {
        Write-Err "pip install completed but 'vastai' is still not on PATH."
        if ($userScripts) {
            Write-Err "Expected at: $(Join-Path $userScripts 'vastai.exe')"
            Write-Err "Test-Path result: $(Test-Path (Join-Path $userScripts 'vastai.exe'))"
        }
        exit 1
    }
}
Write-OK "vastai at $($vastCmd.Source)"

# 2. Vast.ai API key
Write-Step "Configuring Vast.ai API key..."
$existingKeyFile = Join-Path $env:USERPROFILE ".config\vastai\vast_api_key"
if (-not $VastApiKey) {
    if (Test-Path $existingKeyFile) {
        Write-OK "already set in $existingKeyFile (skipping)"
        $env:VAST_API_KEY = (Get-Content $existingKeyFile -Raw).Trim()
    } else {
        $sec = Read-Host "Paste your Vast.ai API key (from https://cloud.vast.ai/cli/)" -AsSecureString
        $bstr = [System.Runtime.InteropServices.Marshal]::SecureStringToBSTR($sec)
        $VastApiKey = [System.Runtime.InteropServices.Marshal]::PtrToStringAuto($bstr)
        [System.Runtime.InteropServices.Marshal]::ZeroFreeBSTR($bstr)
    }
}
if ($VastApiKey) {
    vastai set api-key $VastApiKey | Out-Null
    $env:VAST_API_KEY = $VastApiKey
    Write-OK "VAST_API_KEY exported into this shell"
}

# 3. Dashboard API key
Write-Step "Configuring CASCADE_API_KEY (dashboard auth)..."
$env:CASCADE_API_KEY = $CascadeApiKey
if ($CascadeApiKey -eq "disabled") {
    Write-Warn "set to 'disabled' -- fine for localhost dev; DO NOT use in production"
} else {
    Write-OK "set to a real key"
}
Write-Warn "Browser UI must match: open DevTools and run localStorage.setItem('cascade_api_key', '$CascadeApiKey')"

# 4. SSH key
Write-Step "Configuring SSH key for Vast.ai..."
$candidates = @(
    (Join-Path $env:USERPROFILE ".ssh\id_ed25519"),
    (Join-Path $env:USERPROFILE ".ssh\id_rsa")
)
if (-not $SshKeyPath) {
    foreach ($c in $candidates) { if (Test-Path $c) { $SshKeyPath = $c; break } }
}
if ($SshKeyPath -and (Test-Path $SshKeyPath)) {
    $env:VAST_SSH_KEY_PATH = $SshKeyPath
    Write-OK "VAST_SSH_KEY_PATH=$SshKeyPath"
    $pub = "$SshKeyPath.pub"
    if (Test-Path $pub) {
        Write-Warn "Copy this PUBKEY into https://cloud.vast.ai/account/ -> SSH keys (one-time):"
        Get-Content $pub | ForEach-Object { Write-Host "        $_" -ForegroundColor DarkGray }
    }
} else {
    Write-Warn "no SSH key found at $($candidates -join ' or ')"
    $make = Read-Host "Create one now with ssh-keygen? [y/N]"
    if ($make -eq "y") {
        $target = Join-Path $env:USERPROFILE ".ssh\id_ed25519"
        $sshDir = Split-Path $target -Parent
        if (-not (Test-Path $sshDir)) { New-Item -ItemType Directory -Path $sshDir -Force | Out-Null }
        ssh-keygen -t ed25519 -f $target -N '""' -C "cascade-vastai"
        $env:VAST_SSH_KEY_PATH = $target
        Write-OK "Created $target -- add the .pub to your Vast.ai account."
        $pub = "$target.pub"
        if (Test-Path $pub) {
            Write-Warn "Copy this PUBKEY into https://cloud.vast.ai/account/ -> SSH keys:"
            Get-Content $pub | ForEach-Object { Write-Host "        $_" -ForegroundColor DarkGray }
        }
    } else {
        Write-Warn "Skipping -- you won't be able to launch runs until VAST_SSH_KEY_PATH is set."
    }
}

# 5. Smoke test
Write-Step "Smoke-testing 'vastai search offers'..."
$out = & vastai search offers "gpu_name=A100_SXM4 num_gpus=1 rentable=true verified=true" --on-demand --raw --limit 3 -o dph_total 2>&1
if ($LASTEXITCODE -eq 0) {
    try {
        $parsed = $out | ConvertFrom-Json
        $n = @($parsed).Count
        Write-OK "OK -- Vast.ai returned $n offers. Credentials work."
    } catch {
        Write-Warn "vastai returned exit 0 but output wasn't JSON. Raw output:"
        Write-Host $out -ForegroundColor DarkGray
    }
} else {
    Write-Err "FAILED -- vastai exited $LASTEXITCODE. Output:"
    Write-Host $out -ForegroundColor DarkGray
}

Write-Host ""
Write-Host "Done. Now launch the backend in THIS SAME shell so it inherits the env:" -ForegroundColor Cyan
Write-Host "    uvicorn dashboard_backend.main:app --reload" -ForegroundColor White
Write-Host "(or whatever your dev launcher is -- these env vars are set for the current process only)" -ForegroundColor DarkGray
