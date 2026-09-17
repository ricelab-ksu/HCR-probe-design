#!/bin/bash
# ============================================================================
# bundle_r_mac.sh — build a fully portable, self-contained R runtime for the
# HCR Probe Designer Electron app.
#
# Steps:
#   1. Download the official CRAN R .pkg for the requested macOS architecture
#      (x86_64 or arm64) — ~170 MB.
#   2. Extract the installer payload with `pkgutil --expand-full` (no admin).
#   3. Move the R.framework into resources/r/ next to the app payload.
#   4. Relocate it with fix_r_framework_mac.sh so it runs from the bundle
#      without /Library/Frameworks (this is the cross-platform "trick").
#   5. Install the app's R package dependencies (shiny + deps) directly into
#      the relocated framework's library from CRAN macOS binary packages
#      (no compilation needed).
#
# Usage:
#   bundle_r_mac.sh [version] [arch]
#     version  default: 4.6.0
#     arch     default: detected architecture (arm64 or x86_64)
#
# Tested with R 4.6.0 on macOS arm64 + x86_64 (big-sur base binaries).
# ============================================================================
set -euo pipefail

VERSION="${1:-4.6.0}"
ARCH="${2:-$(uname -m)}"            # arm64 | x86_64
case "$ARCH" in
  arm64|aarch64) ARCH=arm64; SUB=binary ;;   # big-sur-arm64
  x86_64|amd64)  ARCH=x86_64; SUB=binary ;;  # big-sur-x86_64
  *) echo "Unsupported arch: $ARCH" >&2; exit 1 ;;
esac

DIST=$(cd "$(dirname "$0")/.." && pwd)
RES_DIR="$DIST/resources/r"
BUILD_DIR="$DIST/build"

PKG="R-${VERSION}-${ARCH}.pkg"
URL="https://cloud.r-project.org/bin/macosx/big-sur-${ARCH}/base/${PKG}"

echo "==> HCR Probe Designer — R bundle (v${VERSION}, ${ARCH})"
mkdir -p "$BUILD_DIR" "$RES_DIR"

# ---------------------------------------------------------------------------
# 1. Download
# ---------------------------------------------------------------------------
if [ ! -f "$BUILD_DIR/$PKG" ]; then
  echo "==> Downloading $URL"
  curl -fL --retry 3 --progress-bar -o "$BUILD_DIR/$PKG" "$URL"
else
  echo "==> Reusing cached $PKG"
fi
PKG_SIZE=$(du -sh "$BUILD_DIR/$PKG" | awk '{print $1}')
echo "    ($PKG_SIZE)"

# ---------------------------------------------------------------------------
# 2. Extract installer payload (pkgutil works without sudo)
# ---------------------------------------------------------------------------
EXTRACT="$BUILD_DIR/expanded-$ARCH"
if [ ! -e "$EXTRACT/R-fw.pkg/Payload" ] && [ ! -e "$EXTRACT/R.pkg/Payload" ]; then
  rm -rf "$EXTRACT"
  pkgutil --expand-full "$BUILD_DIR/$PKG" "$EXTRACT"
fi

SRC_FW=""
for cand in \
  "$EXTRACT/R-fw.pkg/Payload/R.framework" \
  "$EXTRACT/R.pkg/Payload/Library/Frameworks/R.framework" \
  "$EXTRACT/R.pkg/Payload/R.framework"; do
  if [ -d "$cand" ]; then SRC_FW="$cand"; break; fi
done
if [ -z "$SRC_FW" ]; then
  SRC_FW=$(find "$EXTRACT" -maxdepth 6 -name 'R.framework' -type d 2>/dev/null | head -1)
fi
if [ -z "${SRC_FW:-}" ] || [ ! -d "$SRC_FW" ]; then
  echo "ERROR: R.framework not found in expanded payload" >&2
  find "$EXTRACT" -maxdepth 5 -name 'R.framework' 2>/dev/null || true
  exit 1
fi
echo "==> Found R.framework at $SRC_FW"

# ---------------------------------------------------------------------------
# 3. Stage into resources/r (skip if a complete staging already exists)
# ---------------------------------------------------------------------------
DEST_FW="$RES_DIR/R.framework"
EXISTING_RES=$(ls -d "$DEST_FW"/Versions/*/Resources 2>/dev/null | head -1)
if [ -n "$EXISTING_RES" ] && [ -x "$EXISTING_RES/bin/exec/R" ]; then
  echo "==> Reusing already-staged R.framework at $DEST_FW"
else
  rm -rf "$DEST_FW"
  mkdir -p "$RES_DIR"
  cp -R "$SRC_FW" "$DEST_FW"
  echo "==> Staged R.framework at $DEST_FW"
fi

# Drop AppleDouble residue (._* files) that external volumes sometimes add —
# they break copy steps during packaging.
find "$RES_DIR" -name '._*' -delete 2>/dev/null || true
find "$DEST_FW" -name '._*' -delete 2>/dev/null || true

# ---------------------------------------------------------------------------
# 4. Relocate + re-sign
# ---------------------------------------------------------------------------
bash "$DIST/scripts/fix_r_framework_mac.sh" "$DEST_FW"

# ---------------------------------------------------------------------------
# 5. Install application packages into the bundled library.
#    CRAN macOS "binary" .tgz packages mean shiny + deps land without a
#    compiler. We use the relocated exec/R with R_HOME set explicitly.
# ---------------------------------------------------------------------------
RES="$DEST_FW/Versions/$VERSION/Resources"
[ -d "$RES" ] || RES=$(ls -d "$DEST_FW"/Versions/*/Resources 2>/dev/null | head -1)
execR="$RES/bin/exec/R"
if [ ! -x "$execR" ]; then
  # R 4.6.0+ layout keeps exec/R under bin/exec; fall back to bin/R executable
  execR="$RES/bin/R"
fi
if [ ! -x "$execR" ]; then
  echo "ERROR: cannot find R executable under $RES/bin" >&2
  exit 1
fi

cat > "$BUILD_DIR/install_app_pkgs.R" <<RPKG
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("shiny", "httpuv", "bslib", "sass", "htmltools", "jsonlite")
new_lib <- Sys.getenv("R_LIBS", unset = ".")
dir.create(new_lib, recursive = TRUE, showWarnings = FALSE)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
cat("Installing", length(missing), "missing packages into", new_lib, "\n")
cat(paste(missing, collapse = ", "), "\n")
if (length(missing)) {
  install.packages(missing, lib = new_lib, type = "binary")
}
cat("Session info:\n")
print(sessionInfo())
RPKG

echo "==> Installing R packages (shiny + deps) into bundled library"
R_HOME="$RES" \
R_LIBS="$RES/library" \
R_LIBS_USER="$RES/library" \
R_SHARE_DIR="$RES/share" \
R_INCLUDE_DIR="$RES/include" \
R_DOC_DIR="$RES/doc" \
  "$execR" --vanilla -f "$BUILD_DIR/install_app_pkgs.R"

echo ""
echo "==> Bundle complete."
echo "    R runtime : $RES_DIR/R.framework"
echo "    R libs    : $RES/library"
echo ""
echo "Next: npm install && npm run dist:mac (or npm run dist:mac -- --arm64/--x64)"

# ---------------------------------------------------------------------------
# Optional: verify the relocated framework independently of the system R.
# ---------------------------------------------------------------------------
echo "==> Self-test (port handshake simulator)"
SMOKE="$BUILD_DIR/smoke.R"
cat > "$SMOKE" <<'RSM'
cat("R_HOME=", R.home(), "\n", sep="")
cat("libPaths:", paste(.libPaths(), collapse="; "), "\n")
stopifnot(requireNamespace("shiny", quietly = TRUE))
cat(paste0("HCR_PORT=", httpuv::randomPort(30000L, 60000L)), "\n")
RSM
R_HOME="$RES" R_LIBS="$RES/library" R_LIBS_USER="$RES/library" \
  "$execR" --vanilla -f "$SMOKE" 2>&1 | sed 's/^/    | /'
echo "==> Self-test finished (relocated framework runs standalone)."