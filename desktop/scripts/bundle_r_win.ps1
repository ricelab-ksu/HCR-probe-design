# ============================================================================
# bundle_r_win.ps1 — build a portable R runtime for the Windows build of the
# HCR Probe Designer app.
#
# Rationale: R for Windows is relocatable by design (Rscript.exe computes its
# own R_HOME from its location — that is how R-Portable works), so bundling the
# framework is just: extract the official installer into resources\r, install
# the app's packages into the bundled library, and set R_HOME via env.
#
# Steps:
#   1. Download R-x.y-win.exe from CRAN (Inno Setup based installer).
#   2. Install silently into resources\r (no admin required with user /DIR).
#   3. Install shiny + deps into the bundled library from CRAN binary .zip
#      packages (no compilation).
#   4. Self-test the bundled Rscript.
#
# Usage (PowerShell, from the project root):
#   powershell -ExecutionPolicy Bypass -File scripts\bundle_r_win.ps1 [-Version 4.6.1]
# ============================================================================
param(
    [string]$Version = "4.6.1"
)

$ErrorActionPreference = "Stop"
$Dist = Split-Path -Parent $PSScriptRoot
$ResDir = Join-Path $Dist "resources\r"
$BuildDir = Join-Path $Dist "build"
New-Item -ItemType Directory -Force -Path $BuildDir, $ResDir | Out-Null

$Installer = "R-$Version-win.exe"
# CRAN serves the current Windows R installer directly under base/ (no per-
# version subdir). e.g. https://cloud.r-project.org/bin/windows/base/R-4.6.1-win.exe
$Url = "https://cloud.r-project.org/bin/windows/base/$Installer"
$InstallerPath = Join-Path $BuildDir $Installer

Write-Host "==> HCR Probe Designer - portable R bundle (v$Version, Windows x64)"
if (-not (Test-Path $InstallerPath)) {
    Write-Host "==> Downloading $Url"
    Invoke-WebRequest -Uri $Url -OutFile $InstallerPath -UseBasicParsing
} else {
    Write-Host "==> Reusing cached $Installer"
}

# Silent installer run. /MERGETASKS= avoids desktop icons etc.
Write-Host "==> Extracting R silently into $ResDir"
$targetDir = Join-Path $ResDir "R-$Version"
if (-not (Test-Path (Join-Path $targetDir "bin\Rscript.exe"))) {
    Start-Process -FilePath $InstallerPath -ArgumentList @(
        "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART", "/SP-",
        "/DIR=`"$targetDir`"", "/MERGETASKS=mainr"
    ) -Wait
}
if (-not (Test-Path (Join-Path $targetDir "bin\Rscript.exe"))) {
    throw "R installer did not produce bin\Rscript.exe under $targetDir"
}

$Rscript = Join-Path $targetDir "bin\Rscript.exe"
$RLibDest = Join-Path $ResDir "library"
$RLibUnix = $RLibDest.Replace('\', '/')
New-Item -ItemType Directory -Force -Path $RLibDest | Out-Null

# Packages to install from CRAN binary repo
$installScript = @"
options(repos = c(CRAN = "https://cloud.r-project.org"))
.libPaths(c("$RLibUnix", .libPaths()))
pkgs <- c("shiny", "httpuv", "bslib", "sass", "htmltools", "jsonlite")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
cat("Installing ", length(missing), " missing packages...\n", sep = "")
if (length(missing)) {
  install.packages(missing, lib = "$RLibUnix", type = "binary")
}
print(sessionInfo())
"@
$scriptPath = Join-Path $BuildDir "install_app_pkgs_win.R"
Set-Content -Path $scriptPath -Value $installScript -Encoding ASCII
Write-Host "==> Installing shiny + deps into bundled library"
& $Rscript --vanilla $scriptPath

# Self-test
Write-Host "==> Self-test"
$smoke = Join-Path $BuildDir "smoke_win.R"
Set-Content -Path $smoke -Value @"
.libPaths(c("$RLibUnix", .libPaths()))
cat("R_HOME=", R.home(), "\n", sep="")
stopifnot(requireNamespace("shiny", quietly = TRUE))
cat(paste0("HCR_PORT=", httpuv::randomPort(30000L, 60000L)), "\n")
"@ -Encoding ASCII
$env:R_LIBS = $RLibDest
& $Rscript --vanilla $smoke

Write-Host ""
Write-Host "==> Windows R bundle complete."
Write-Host "    R runtime : $targetDir"
Write-Host "    R libs    : $RLibDest"
Write-Host ""
Write-Host "Next: npm install && npm run dist:win"